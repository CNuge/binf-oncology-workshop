rule sync_consent_information:
    """Sync the jsons that contain consent status locally

    Note that we've flagged the output as temp here.
    That means that every time the pipeline is run,
    the consent jsons will be pulled fresh and the consent
    files re-generated.

    Also note the choice to use aws s3 sync instead of
    snakemake s3 remotes - s3 remotes become non-performant
    above a couple dozen files.

    NB each individual has their own json object on the s3 bucket.
    """
    output:
        temp("results/sync/consent_jsons.done.flag"),
    params:
        bucket=config["s3_input_targets"]["consent_bucket"],
    shell:
        "mkdir -p input/{params.bucket} && "
        'aws s3 sync --exclude "*" --include "*.json" '
        "{params.bucket}/ input/{params.bucket}/ "
        "--only-show-errors &> {output}"


rule list_consent_jsons:
    input:
        flag="results/sync/consent_jsons.done.flag",
    output:
        "results/consent_lists/file_list.tsv",
    params:
        json_dir="input/{}".format(config["s3_input_targets"]["consent_bucket"]),
    shell:
        "for i in $(ls {params.json_dir}); do echo {params.json_dir}/$i >> {output}; done"


rule generate_consent_lists:
    """Generates two tsv files in consent_lists folder. 

    This rule calls the script parse_consent.py which parses s3 json into tidy pandas dataframes.

    outputs are:
    'consented.tsv' file will contain the individuals that are okay to use, 
    'unconsented.tsv' are individuals that should be filtered out. 
    """
    input:
        filelist="results/consent_lists/file_list.tsv",
        flag="results/sync/consent_jsons.done.flag",
    output:
        "results/consent_lists/consented.tsv",
        "results/consent_lists/unconsented.tsv",
    params:
        study_name=config["study_name_for_consent"],
        out_prefix="results/consent_lists",
        json_dir="input/{}".format(config["s3_input_targets"]["consent_bucket"]),
    shell:
        "python workflow/scripts/parse_consent.py {input.filelist} '{params.study_name}' {params.out_prefix} && "
        "rm -r {params.json_dir}"


rule apply_consent:
    """
    Match individuals in the manifest to their consent status.

    Produces a new manifest with unconsented individuals removed,
    should any be found in the data.

    Uses the seqbiopy apply_consent CLI tool.
    Note that seqbiopy also contains a parse_consent CLI tool (to replace preceeding rule)
    but the use of both a script file and seqbiopy tool are shown in the cookie cutter
    for demonstration purposes. 
    """
    input:
        # NOTE - these required inputs ensure family structure info
        # and raw consent info are acquired and parsed beforehand
        consented_file="results/consent_lists/consented.tsv",
        unconsented_file="results/consent_lists/unconsented.tsv",
        manifest_file=config["manifest"],
    output:
        consented_mainfest="results/consent_lists/consented_mainfest.tsv",
        manifest_check_flag="results/flags/consented_check.flag",
    params:
        manifest_id="uuid",
        consent_id="global_participant_id",
    shell:
        " apply_consent -m {input.manifest_file} -c {input.consented_file} "
        " -f {output.manifest_check_flag} -o {output.consented_mainfest} "
        " -i {params.manifest_id} -n {params.consent_id} "
