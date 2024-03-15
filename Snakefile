METHODS_QQ = ["combat","harmony","mnn","seurat","single"]
METHODS = ["","bbknn","combat","harmony","mnn","scvi","seurat","seurat_v2"]
ORIG_OUTPUT = [".h5ad","_single.h5ad","_to_csv.csv","_to_csv_raw.csv","_variable_genes.csv","_pca_to_csv.csv","_batch.csv"]
data_dir={
    'pbmc4k':"data/filtered_gene_bc_matrices/hg19/",
    'neuro': "data/1M_neurons_neuron20k.h5",


}
heart_files=[f"data/heart-{i}.txt" for i in range(25)]
neuro_files=[f"data/neuro-{i}.txt" for i in range(25)]
pbmc3k_files=[f"data/pbmc3k-{i}.txt" for i in range(25)]
pbmc4k_files=[f"data/pbmc4k-{i}.txt" for i in range(25)]
diffexp_simul_pbmc_files=[f"data/diffexp-simul-pbmc-{i}.txt" for i in range(25)]
diffexp_simul_neuro_files=[f"data/diffexp-simul-neuro-{i}.txt" for i in range(25)]
diffexp_pbmc_files=[f"data/diffexp-pbmc-{i}.txt" for i in range(25)]
diffexp_neuro_files=[f"data/diffexp-neuro-{i}.txt" for i in range(25)]
simul_pbmc_files =[f"data/simul-pbmc-{i}.txt" for i in range(25)]
simul_neuro_files =[f"data/simul-neuro-{i}.txt" for i in range(25)]
all_files = heart_files + neuro_files + pbmc3k_files + pbmc4k_files + diffexp_simul_pbmc_files + diffexp_simul_neuro_files + diffexp_pbmc_files + diffexp_neuro_files + simul_pbmc_files + simul_neuro_files


def make_output_names(prefix,simul=False,simul_string=None,prepend_letters=False):
    postfixes = ".h5ad _single.h5ad _to_csv.csv _to_csv_raw.csv _variable_genes.csv _pca_to_csv.csv _batch.csv".split()
    if simul:
        if simul_string is not None:
            ss = "simul-"+simul_string
        else:
            ss = "simul-"+prefix
        ls = [
            rf"data/{{adata, {prefix}(?![\w\d_-])|{prefix}-\d+|{ss}-\d+}}{pf}"
            for pf in postfixes
        ]
        #ls[-1] = ls[-1][:-1]
    else:
        ls = [
            rf"data/{{adata, {prefix}(?![\w\d_-])|{prefix}-\d+}}{pf}"
            for pf in postfixes
        ]
        ls[-1] = ls[-1][:-1]
    if prepend_letters:
        letters = ["a","b","c","d","e","f","g","h"]
        for i in range(len(ls)):
            ls[i] = letters[i] + " = " + ls[i]
    return ls


 
rule all:
    input:
        all_files
    output:
        "data/all.txt"
    shell:
        "echo job done > {output}"

rule run_all_pbmc3k:
    input: pbmc3k_files
    output:
        "data/pbmc3k-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_neuro:
    input: neuro_files
    output:
        "data/neuro-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_hearts:
    input: heart_files
    output:
        "data/heart-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_pbmc4k:
    input: pbmc4k_files
    output:
        "data/pbmc4k-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_diffexp_pbmc:
    input: diffexp_pbmc_files
    output:
        "data/diffexp_pbmc-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_diffexp_neuro:
    input: diffexp_neuro_files
    output:
        "data/diffexp-neuro-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_diffexp_simul_pbmc:
    input: diffexp_simul_pbmc_files
    output:
        "data/diffexp-simul-pbmc-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_diffexp_simul_neuro:
    input: diffexp_simul_neuro_files
    output:
        "data/diffexp-simul-neuro-all.txt"
    shell:
        "echo job done > {output}"
rule run_all_simul_pbmc:
    input: simul_pbmc_files
    output:
        "data/simul-pbmc-all.txt"
rule run_all_simul_neuro:
    input: simul_neuro_files
    output:
        "data/simul-neuro-all.txt"
    shell:
        "echo job done > {output}"




rule create_adata_pbmc_orig:
    input:
        "data/filtered_gene_bc_matrices/hg19/"
    output:
        make_output_names("pbmc3k")
    params:
        input="{adata}"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_base.py"

#expand("data/{{adata,adata-simul-pbmc(?![\w\d_-])|adata-simul-pbmc-\d+}}{out}",out=ORIG_OUTPUT)
rule create_simul_pbmc:
    input:
        "data/filtered_gene_bc_matrices/hg19/"
    output:
        make_output_names("diffexp-simul-pbmc",True,"pbmc")

        
    params:
        input="{adata}"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_simulated.py"


rule create_simul_neuro:
    input:
        "data/1M_neurons_neuron20k.h5"
    output:
        make_output_names("diffexp-simul-neuro",True,"neuro")
    params:
        input="{adata}"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_simulated.py"


rule create_adata_neuro:
    input:
        "data/1M_neurons_neuron20k.h5"
    output:
        make_output_names("neuro")

    conda:
        "envs/env.yml"
    script:
        "scripts/mod_base_neuro.py"


#DIFFEXP H5AD creation

rule create_adata_pbmc_diffexp:
    input:
        "data/filtered_gene_bc_matrices/hg19/"
    output:
        make_output_names("diffexp-pbmc")
    params:
        input="{adata}"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_base.py"

rule create_adata_neuro_diffexp:
    input:
        "data/1M_neurons_neuron20k.h5"
    output:
        make_output_names("diffexp-neuro")

    conda:
        "envs/env.yml"
    script:
        "scripts/mod_base_neuro.py"









#------------------------ Resampling






########## SIMULATED AS METHOD


rule create_adata_heart:
    input:
        "data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R1_001.fastq.gz",
        "data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R2_001.fastq.gz",
        "data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R1_001.fastq.gz",
        "data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R2_001.fastq.gz"

    output:
        a = "data/{adata,heart(?![\w\d_-])|heart-\d+}.h5ad",
        b = "data/{adata,heart(?![\w\d_-])|heart-\d+}_single.h5ad",
        c = "data/{adata,heart(?![\w\d_-])|heart-\d+}_to_csv.csv",
        d = "data/{adata,heart(?![\w\d_-])|heart-\d+}_to_csv_raw.csv",
        e = "data/{adata,heart(?![\w\d_-])|heart-\d+}_variable_genes.csv",
        f = "data/{adata,heart(?![\w\d_-])|heart-\d+}_pca_to_csv.csv",
        g = "data/{adata,heart(?![\w\d_-])|heart-\d+}_batch.csv"
        #make_output_names("heart",prepend_letters=True)
    params:
        input="{adata}"
    conda:
        "envs/env_bus.yml"
    script:
        "scripts/mod_base_heart.py"



rule create_adata_pbmc4k:
    input:
        "data/fastqs/pbmc4k_S1_L001_R1_001.fastq.gz",
        "data/fastqs/pbmc4k_S1_L001_R2_001.fastq.gz",
        "data/fastqs/pbmc4k_S1_L002_R1_001.fastq.gz",
        "data/fastqs/pbmc4k_S1_L002_R2_001.fastq.gz"
    output:
        a = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}.h5ad",
        b = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_single.h5ad",
        c = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_to_csv.csv",
        d = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_to_csv_raw.csv",
        e = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_variable_genes.csv",
        f = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_pca_to_csv.csv",
        g = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_batch.csv"
        #make_output_names("pbmc4k",prepend_letters=True)
    params:
        input="{adata}"
    conda:
        "envs/env_bus.yml"
    script:
        "scripts/mod_base_pbmc4k.py"


rule pbmc4k_downsample_run:
    input:
        rules.create_adata_pbmc4k.output.a
    output:
        a = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_down.h5ad"
    priority: 1
    conda:
        "envs/env_down.yml" 
    script:
        "scripts/mod_down.py"


rule pbmc4k_resample_run:
    input:
        rules.create_adata_pbmc4k.output.a
    output:
        a = "data/{adata,pbmc4k(?![\w\d_-])|pbmc4k-\d+}_resample.h5ad"
    conda:
        "envs/env_bus.yml" 
    script:
        "scripts/mod_resample.py"



rule heart_downsample_run:
    input:
        rules.create_adata_heart.output.a
    output:
        a = "data/{adata,heart(?![\w\d_-])|heart-\d+}_down.h5ad"
    priority: 1
    conda:
        "envs/env_down.yml" 
    script:
        "scripts/mod_down.py"


rule heart_resample_run:
    input:
        rules.create_adata_heart.output.a
    output:
        a = "data/{adata,heart(?![\w\d_-])|heart-\d+}_resample.h5ad"
    conda:
        "envs/env_bus.yml" 
    script:
        "scripts/mod_resample.py"

















rule mnn_run_1:
    input:
        "data/{adata}_to_csv.csv",
        "data/{adata}_batch.csv",

    output:
        temp("data/r_to_{adata}_mnn.csv")
    conda:
        "envs/env_mnn.yml" 
    script:
        "scripts/mod_mnn.R"
rule mnn_run_2:
    input:
        "data/{adata}.h5ad",
        "data/r_to_{adata}_mnn.csv"
    output:
        "data/{adata}_mnn.h5ad",
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_mnn.py"


rule seurat_run_1:
    input:
        "data/{adata}_to_csv.csv",
        "data/{adata}_batch.csv",
        "data/{adata}_variable_genes.csv"

        #"data/variable_genes.csv"
    output:
        temp("data/r_to_{adata}_seurat.csv")
    params:
        mode = "FALSE"
    conda:
        "envs/env_seurat.yml"
    script:
        "scripts/mod_seurat.R"
rule seurat_run_2:
    input:
        "data/{adata}.h5ad",
        "data/r_to_{adata}_seurat.csv"
    output:
        "data/{adata}_seurat.h5ad",
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_seurat.py"

rule seuratv2_run_1:
    input:
        "data/{adata}_to_csv_raw.csv",
        "data/{adata}_batch.csv",
        "data/{adata}_variable_genes.csv"
    output:
        temp("data/r_to_{adata}_seurat_v2.csv")
    conda:
        "envs/env_seurat.yml"
    params:
        mode = "TRUE"
    script:
        "scripts/mod_seurat.R"
rule seuratv2_run_2:
    input:
        "data/{adata}.h5ad",
        "data/r_to_{adata}_seurat_v2.csv"
    output:
        "data/{adata}_seurat_v2.h5ad",
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_seurat_v2.py"

rule harmony_run_1:
    input:
        "data/{adata}_to_csv.csv",
        "data/{adata}_pca_to_csv.csv",
        "data/{adata}_batch.csv"
    output:
        temp("data/r_to_{adata}_harmony.csv")
    conda:
        "envs/env_harmony.yml"
    script:
        "scripts/mod_harmony.R"
rule harmony_run_2:
    input:
        "data/{adata}.h5ad",
        "data/r_to_{adata}_harmony.csv"
    output:
        "data/{adata}_harmony.h5ad"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_harmony.py"

rule combat_run:
    input:
        "data/{adata}.h5ad"
    output:
        "data/{adata}_combat.h5ad"
    conda:
        "envs/env.yml"
    script:
        "scripts/mod_combat.py"

rule bbknn_run:
    input:
        "data/{adata}.h5ad"
    output:
        "data/{adata}_bbknn.h5ad"
    conda:
        "envs/env_bbknn.yml"
    script:
        "scripts/mod_bbknn.py"

rule scvi_run:
    input:
        "data/{adata}.h5ad"
    output:
        "data/{adata}_scvi.h5ad"
    conda:
        "envs/env_scvi.yml"
    script:
        "scripts/mod_scvi.py"
rule liger_run:
    input:
        "data/{adata}.h5ad"
    output:
        "data/{adata}_liger.h5ad"
    conda:
        "envs/env_liger.yml"
    script:
        "scripts/mod_pyliger.py"
rule liger_v2_run:
    input:
        "data/{adata}.h5ad"
    output:
        "data/{adata}_liger_v2.h5ad"
    conda:
        "envs/env_liger.yml"
    script:
        "scripts/mod_pyliger_v2.py"
    





rule run_pbmc_orig:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_liger.h5ad",
        "data/{adata}_liger_v2.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad"
    output:
        "data/{adata,adata(?![\w\d_-])|adata-\d+}.txt"

    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"


rule run_neuro:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_liger.h5ad",
        "data/{adata}_liger_v2.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad"

    output:
        "data/{adata,neuro|neuro-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"


rule run_pbmc4k:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_liger.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad",
        #"data/{adata}_down.h5ad",
        #"data/{adata}_resample.h5ad"
        rules.pbmc4k_downsample_run.output.a,
        rules.pbmc4k_resample_run.output.a
    output:
        "data/{adata,pbmc4k|pbmc4k-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"




rule run_heart:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_liger.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad",
        #"data/{adata}_down.h5ad",
        rules.heart_downsample_run.output.a,
        rules.heart_resample_run.output.a
        #"data/{adata}_resample.h5ad"
    output:
        "data/{adata,heart|heart-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"


####### SIMULATED


rule run_simul_pbmc:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad"
    output:
        "data/{adata,simul-pbmc|simul-pbmc-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"

rule run_simul_neuro:
    input:
        #expand("{{adata}}_{method}.h5ad",method=METHODS)
        #expand("data/{{adata}}.h5ad", i = ITER),
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad",
        "data/{adata}_seurat_v2.h5ad"
    output:
        "data/{adata,simul-neuro|simul-neuro-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output.py"






#DIFFERENTIAL EXPRESSION
rule run_diffexp_pre_r:
    input:
        #"data/{adata}_to_csv.csv",
        #"data/{adata}_batch.csv",
        "data/{adata}.h5ad",
        "data/{adata}_bbknn.h5ad",
        "data/{adata}_combat.h5ad",
        "data/{adata}_harmony.h5ad",
        "data/{adata}_liger.h5ad",
        "data/{adata}_mnn.h5ad",
        "data/{adata}_scvi.h5ad",
        "data/{adata}_seurat.h5ad"
        #"data/{adata}_seurat_v2.h5ad"
        #"data/{adata}_leiden.txt",
        #"data/{adata}_clusters.csv"
    output:
        a = "data/{adata}_batch_p_vals.csv",
        b = "data/{adata}_orig_p_vals.csv"
    conda:
        "envs/env_mast.yml"
    script:
        "scripts/mod_mast.R"


rule run_diffexp_pbmc:
    input:
        rules.run_diffexp_pre_r.output.a,
        rules.run_diffexp_pre_r.output.b

    output:
        "data/{adata,diffexp-pbmc|diffexp-pbmc-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output_diffexp.py"


rule run_diffexp_neuro:
    input:
        rules.run_diffexp_pre_r.output.a,
        rules.run_diffexp_pre_r.output.b

    output:
        "data/{adata,diffexp-neuro|diffexp-neuro-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output_diffexp.py"


rule run_diffexp_pbmc_simul:
    input:
        rules.run_diffexp_pre_r.output.a,
        rules.run_diffexp_pre_r.output.b

    output:
        "data/{adata,diffexp-simul-pbmc|diffexp-simul-pbmc-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output_diffexp.py"



rule run_diffexp_pbmc_neuro:
    input:
        rules.run_diffexp_pre_r.output.a,
        rules.run_diffexp_pre_r.output.b

    output:
        "data/{adata,diffexp-simul-neuro|diffexp-simul-neuro-\d+}.txt"
    conda:
        "envs/env.yml"
    script:
        "scripts/create_output_diffexp.py"