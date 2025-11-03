#!/usr/bin/env python3

import argparse
import json

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description="Take a list of per-subject jsons and create one rectangular tidy tsv containing global and study-specific subject IDs as well as enrollment (aka consent) status"
    )
    parser.add_argument(
        "filename_list",
        type=str,
        help="list of json files containing per-subject consent info",
    )
    parser.add_argument(
        "study_name",
        type=str,
        help="name of the study for which consent is being evaluated",
    )
    parser.add_argument(
        "out_prefix",
        type=str,
        help="prefix of output tsv file (will be appended with <yaml key>.tsv",
    )
    results = parser.parse_args()
    return (results.filename_list, results.study_name, results.out_prefix)


if __name__ == "__main__":
    filename_list, study_name, outprefix = get_args()
    json_files = []
    with open(filename_list, "r") as j:
        json_files = j.readlines()

    df = pd.DataFrame()
    for filename in json_files:
        with open(filename.rstrip(), "r") as f:
            subject_data = json.load(f)

        for study in subject_data["studies"]:
            if study["name"].lower() == study_name.lower() and "enrollment_status" in study.keys():
                tmp_df = pd.DataFrame(
                    {
                        "global_participant_id": subject_data["participant_id"],
                        "global_short_id": subject_data["participant_short_id"],
                        "study_participant_id": study["participant_id"],
                        "study_short_id": study["participant_short_id"],
                        "enrollment_status": study["enrollment_status"],
                    },
                    index=[0],
                )

        df = pd.concat([df, tmp_df])

    # Note - "partially_withdrawn" means individuals have asked to not be contacted further
    # DOES NOT indicate that they cannot be used for certain types of analyses.
    # Note - enrolled individuals may withdraw consent at a later date and move to the
    # "withdrawn" category.
    confirmed_consent = ["enrolled", "partially_withdrawn"]
    # Note - "pending" or "enrolling" individuals may change to enrolled at a later date.
    unconsented = ["pending", "enrolling", "withdrawn"]
    consented_df = df.loc[df["enrollment_status"].isin(confirmed_consent)].drop_duplicates()
    unconsented_df = df.loc[df["enrollment_status"].isin(unconsented)].drop_duplicates()

    consented_df.to_csv(outprefix + "/consented.tsv", sep="\t", index=False)
    unconsented_df.to_csv(outprefix + "/unconsented.tsv", sep="\t", index=False)
