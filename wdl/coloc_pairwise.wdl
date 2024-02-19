version 1.0

workflow ColocSusie{
    input {
        File map1
        File list1
        String bgenInf1
        String clumpGeno1
        Float sigThresh1

        File map2
        String list2
        String bgenInf2
        String clumpGeno2
        Float sigThresh2

        Int windowKB = 1500
    }

    Array[Array[String]] probs1 = read_tsv(list1)
    Array[Array[String]] probs2 = read_tsv(list2)

    scatter(prob1 in probs1){
        scatter(prob2 in probs2){
            String outPrefix = prob1[0] + "---" + prob2[0]
            call mergeClump {input: summary1=prob1[1], map1=map1, sigThresh1=sigThresh1, subsample1=prob1[2], clumpGeno1=clumpGeno1, summary2=prob2[1], map2=map2, sigThresh2=sigThresh2, subsample2=prob2[2], clumpGeno2=clumpGeno2, windowKB=windowKB, outPrefix=outPrefix}
        }
    }

    call mergeRegion {
        input: hadResults=flatten(mergeClump.hadResults), regions=flatten(mergeClump.region), summaries=flatten(mergeClump.commSuma), samples1=flatten(mergeClump.sample1), samples2=flatten(mergeClump.sample2) 
    }

    Array[Array[String]] regions = read_tsv(mergeRegion.list)
    scatter(region in regions){
        String outPrefix2 = region[0] + "___" + region[1]
        call colocPrep{
            input: commSuma=region[7], region=region[1], outPrefix=outPrefix2, bgenInf1=bgenInf1, bgenInf2=bgenInf2
        }

        String outPrefixLD1 = outPrefix2 + "_LD1"
        call colocLD as colocLD1 {
            input: zfile=colocPrep.zfile, bgenGeno=colocPrep.bgenGeno1, bgenIndex=colocPrep.bgenIndex1, bgenSample=colocPrep.bgenSample1, N=colocPrep.N1, subsample=region[8], outPrefix=outPrefixLD1
        }

        String outPrefixFinemap1 = outPrefix2 + "_finemap1"
        call finemap as finemap1 {
            input: ldz=colocLD1.ldz, suma=colocPrep.suma1, outPrefix=outPrefixFinemap1
        }

        String outPrefixLD2 = outPrefix2 + "_LD2"
        call colocLD as colocLD2 {
            input: zfile=colocPrep.zfile, bgenGeno=colocPrep.bgenGeno2, bgenIndex=colocPrep.bgenIndex2, bgenSample=colocPrep.bgenSample2, N=colocPrep.N2, subsample=region[9], outPrefix=outPrefixLD2
        }

        String outPrefixFinemap2 = outPrefix2 + "_finemap2"
        call finemap as finemap2{
            input: ldz=colocLD2.ldz, suma=colocPrep.suma2, outPrefix=outPrefixFinemap2
        }

        String outPrefixColoc = outPrefix2 + "_coloc"
        call coloc {
            input: finemap1=finemap1.res, finemap2=finemap2.res, outPrefix=outPrefixColoc
        }

        call mergeResults {
            input: trait1=region[0], trait2=region[1], suma=region[7], coloc_res=coloc.res, finemap1_res=finemap1.res, finemap2_res=finemap2.res
        }
    }

    call mergeAll {
        input: results=mergeResults.res
    }

    output {
        File merged = mergeAll.res
        File list = mergeRegion.list
        Array[File] finemapped1 = finemap1.res
        Array[File] finemapped2 = finemap2.res
        Array[File] colocs = coloc.res
    }
}

task mergeAll {
    input {
        Array[File] results
    }

    command <<<
        awk 'FNR>1 || NR==1' ~{sep=' ' results} > colocs.tsv
    >>>

    runtime{
        cpu: 2
        memory: "2 GB"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
        disks: "local-disk 10 SSD"
    }

    output {
        File res = "colocs.tsv"
    }
}


task mergeResults {
    input {
        String trait1
        String trait2
        File suma
        File coloc_res
        File finemap1_res
        File finemap2_res
    }
    command <<<
        mergeResults.R "~{trait1}" "~{trait2}" "~{suma}" "~{coloc_res}" "~{finemap1_res}" "~{finemap2_res}"
    >>>

    runtime{
        cpu: 2
        memory: "4 GB"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
        disks: "local-disk 10 SSD"
    }

    output {
        File res = trait1 + "___" + trait2 + "_combined.tsv"
    }
}

task mergeRegion{
    input {
        Array[Boolean] hadResults
        Array[File] regions
        Array[String] summaries
        Array[String] samples1
        Array[String] samples2
    }

    command <<<
        #cat ~{sep=' ' regions} > regions.list
        hadRes=(~{sep=' ' hadResults})
        regions=(~{sep=' ' regions})
        summaries=(~{sep=' ' summaries})
        samples1=(~{sep=' ' samples1})
        samples2=(~{sep=' ' samples2})

        for idx in "${!hadRes[@]}"; do
            if [[ "${hadRes[$idx]}" == "true" ]]; then
                while read region; do
                    echo -e "${region}\t${summaries[$idx]}\t${samples1[$idx]}\t${samples2[$idx]}" >> regions.tsv
                done < <(cat ${regions[$idx]})
            fi
        done

    >>>

    runtime{
        cpu: 2
        memory: "4 GB"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
        disks: "local-disk 20 SSD"
    }

    output {
        File list = "regions.tsv"
    }
}

task mergeClump{
    input {
        File summary1 
        File map1
        Float sigThresh1
        File subsample1
        String clumpGeno1

        File summary2
        File map2
        Float sigThresh2
        File subsample2
        String clumpGeno2

        Float r2 = 0.3
        Int windowKB

        String outPrefix
    }

    command <<<

        # gs file to local file system
        # input 1: gs bucket
        # input 2: local mount folder
        # output: mountpoint localpath
        function gs2local() {
            local gstemp=$1
            local bucket="${gstemp#*//}"
            local bucket="${bucket%%/*}"

            local fullbucket="gs://$bucket"

            local mountpoint="$2/${bucket}"
            local localpath=${gstemp//$fullbucket/$mountpoint}

            echo  ${localpath} ${mountpoint} ${bucket}
        }

        echo "false" > had_results
        echo "=====Tidying summary 1"
        mungeSuma.R -i ~{summary1} -m ~{map1} -o ~{outPrefix}_m1 -s ~{subsample1}

        echo "=====Tidying summary 2"
        mungeSuma.R -i ~{summary2} -m ~{map2} -o ~{outPrefix}_m2 -s ~{subsample2}
        
        echo "=====Merging the two summary"
        mergeSuma.R ~{outPrefix}_m1 ~{outPrefix}_m2 ~{sigThresh1} ~{sigThresh2} ~{outPrefix}

        echo "====Mount genotype1"
        read geno1 mount1 bucket1 < <(gs2local ~{clumpGeno1} /cromwell_root/gcsfuse)
        mkdir -p $mount1
        gcsfuse --implicit-dirs ${bucket1} ${mount1} 

        echo "=====Clumping summary 1"
        plink --bfile ${geno1} --clump ~{outPrefix}_clump1.txt --extract ~{outPrefix}_clump1.snplist --clump-p1 ~{sigThresh1} --clump-p2 ~{sigThresh1} --clump-r2 ~{r2} --clump-kb ~{windowKB} --out ~{outPrefix}_clump1 --memory 12288  --threads 1
        fusermount -uz ${mount1}

        echo "====Mount genotype2"
        read geno2 mount2 bucket2 < <(gs2local ~{clumpGeno2} /cromwell_root/gcsfuse)
        mkdir -p $mount2
        gcsfuse --implicit-dirs ${bucket2} ${mount2} 

        echo "=====Clumping summary 2"
        plink --bfile ${geno2} --clump ~{outPrefix}_clump2.txt --extract ~{outPrefix}_clump2.snplist --clump-p1 ~{sigThresh2} --clump-p2 ~{sigThresh2} --clump-r2 ~{r2} --clump-kb ~{windowKB} --out ~{outPrefix}_clump2 --memory 12288 --threads 1
        fusermount -uz ${mount2}

        echo "=====Merging the ranges"
        mergeRange.R ~{outPrefix} ~{windowKB}

        echo "=====Test results"
        if [ -f "~{outPrefix}.region.txt" ]; then
            echo "Have region file"
            if [ $(wc -l ~{outPrefix}.region.txt | awk '{print $1}') -ge "1" ]; then
                echo " valid regions"
                cat "~{outPrefix}.region.txt"
                echo "true" > had_results
            fi
        fi

        tar -zcvf ~{outPrefix}_clump.tar.gz ~{outPrefix}.raw.clump

    >>>

    runtime{
        cpu: 2
        memory: "16 GB"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
        disks: "local-disk 25 HDD"
    }

    output {
        File region = outPrefix + ".region.txt"
        File commSuma = outPrefix + ".rds"
        File clump = outPrefix + "_clump.tar.gz"
        Boolean hadResults = read_boolean("had_results")
        File sample1 = subsample1
        File sample2 = subsample2
        #Array[File] suma1 = glob("*_m1.gz")
        #Array[File] suma2 = glob("*_m2.gz")
    }
}

task colocPrep{
    input{
        File commSuma
        String region
        String outPrefix
        File bgenInf1
        File bgenInf2
    }

    command <<<
        colocPrep.R --suma ~{commSuma} --region ~{region} --bgenInf1 ~{bgenInf1} --bgenInf2 ~{bgenInf2} --prefix ~{outPrefix}
    >>>

    output {
        File zfile = outPrefix + ".z"
        String bgenGeno1 = read_string(outPrefix + ".bgenGeno1.txt")
        String bgenIndex1 = read_string(outPrefix + ".bgenIndex1.txt")
        String bgenSample1 = read_string(outPrefix + ".bgenSample1.txt")
        String bgenGeno2 = read_string(outPrefix + ".bgenGeno2.txt")
        String bgenIndex2 = read_string(outPrefix + ".bgenIndex2.txt")
        String bgenSample2 = read_string(outPrefix + ".bgenSample2.txt")
 
        Int N1 = read_int(outPrefix + ".N1.txt")
        Int N2 = read_int(outPrefix + ".N2.txt")

        Int nsnps = read_int(outPrefix + ".n.txt")

        File suma1 = outPrefix + ".suma1.rds"
        File suma2 = outPrefix + ".suma2.rds"
    }
    runtime {
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 20 SSD"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
    }
}



task colocLD{
    input{
        String bgenGeno
        File bgenIndex
        File bgenSample
        Int N
        File zfile
        File subsample
        String outPrefix

        Int cpu = 4

        String master = outPrefix + ".master"
    }
    
    command <<<
        # gs file to local file system
        # input 1: gs bucket
        # input 2: local mount folder
        # output: mountpoint localpath
        function gs2local() {
            local gstemp=$1
            local bucket="${gstemp#*//}"
            local bucket="${bucket%%/*}"

            local fullbucket="gs://$bucket"

            local mountpoint="$2/${bucket}"
            local localpath=${gstemp//$fullbucket/$mountpoint}

            echo  ${localpath} ${mountpoint} ${bucket}
        }

        echo "====Mount bgen"
        read bgen mountp bucket < <(gs2local ~{bgenGeno} /cromwell_root/gcsfuse)
        mkdir -p $mountp
        gcsfuse --implicit-dirs ${bucket} ${mountp} 

        awk '
        BEGIN {
            OFS = ";"
            print "z", "bgen", "bgi", "bdose", "bcor", "ld", "sample", "incl", "n_samples"
            print "~{zfile}", "'${bgen}'", "~{bgenIndex}", "~{outPrefix}.bdose", "~{outPrefix}.bcor", "~{outPrefix}.ld", "~{bgenSample}", "~{subsample}", "~{N}"
        }' > ~{master}

        ldstore --in-files ~{master} --write-bcor --write-bdose --bdose-version 1.1 --memory 10 --n-threads ~{cpu} 
        ldstore --in-files ~{master} --bcor-to-text
        bgzip -@ ~{cpu} ~{outPrefix}.ld
        
        echo "==== FILE sizes"
        ls -alhtr

        fusermount -uz $mountp
    >>>

    output{
        File bcor = outPrefix + ".bcor"
        File ldz = outPrefix + ".ld.gz"
    }

    runtime {
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        cpu: "${cpu}"
        memory: "12 GB"
        disks: "local-disk 50 SSD"
        zones: "europe-west1-b"
    }
}

task finemap {
    input{
        File ldz
        File suma
        String outPrefix
        Int mem = 20
    }

    command <<<
        finemap.R --suma ~{suma} --ldz ~{ldz} --output ~{outPrefix}
    >>>

    output{
        File res = outPrefix + ".rds"
    }

    runtime{
        cpu: 4
        memory: "${mem} GB"
        disks: "local-disk 15 SSD"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
    }
}

task coloc {
    input {
        File finemap1
        File finemap2
        String outPrefix
    }

    command <<<
        coloc.R --finemap1 ~{finemap1} --finemap2 ~{finemap2} --output ~{outPrefix}
    >>>

    output {
        File res = outPrefix + ".rds"
    }
    runtime{
        cpu: 1
        memory: "5 GB"
        disks: "local-disk 15 SSD"
        docker: "europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/coloc.susie:0.1.4"
        zones: "europe-west1-b"
    }
}

