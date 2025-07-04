#!/usr/bin/env python3

# based on work by Damien Farrell https://dmnfarrell.github.io/bioinformatics/bcftools-csq-gff-format
import sys

def GFF_bcftools_format(in_handle, out_handle):
    """Convert a bacterial genbank file from NCBI to a GFF3 format that can be used in bcftools csq.
    see https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt#L1066-L1098.
    Args:
        in_file: genbank file
        out_file: name of GFF file
    """

    from BCBio import GFF

    in_handle = open(in_handle)
    out_handle = open(out_handle, "w")
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from copy import copy, deepcopy
    from Bio import SeqIO

    # handle cases with no 'gene' or 'locus_tag qualifiers
    hypothetical_count = 0

    for record in SeqIO.parse(in_handle, "genbank"):
        # make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        # loop over all features
        for feat in record.features:
            q = feat.qualifiers
            # remove some unecessary qualifiers
            for label in ["note", "translation", "product", "experiment"]:
                if label in q:
                    del q[label]

            if feat.type == "CDS":
                print(feat)
                # use the CDS feature to create the new lines
                try:
                    tag = q["gene"][0]  # q['locus_tag'][0]
                except:
                    try:
                        tag = q['locus_tag'][0]
                    except:
                        hypothetical_count += 1
                        tag = f"hypothetical_{hypothetical_count}"
                try:
                    protein_id = q["protein_id"][0]
                except:
                    protein_id = "Unknown"
                q["ID"] = "CDS:%s" % protein_id
                q["biotype"] = "protein_coding"

                for i, new_loc in enumerate(
                    feat.location.parts
                    if (hasattr(feat.location, "parts"))
                    else (feat.location,)
                ):
                    new_feat = deepcopy(feat)
                    tr_id = "transcript:%s" % (protein_id + "_" + str(i))
                    new_feat.qualifiers["Parent"] = tr_id
                    new_feat.location = new_loc
                    new.features.append(new_feat)
                    # create mRNA feature
                    m = SeqFeature(feat.location, type="mRNA")
                    q2 = m.qualifiers
                    q2["ID"] = tr_id
                    q2["Parent"] = "gene:%s" % tag
                    q2["biotype"] = "protein_coding"
                    new.features.append(m)
            elif feat.type == "gene":
                tag = q["gene"][0]
                # edit the gene feature
                q = feat.qualifiers
                q["ID"] = "gene:%s" % tag
                q["biotype"] = "protein_coding"
                q["Name"] = q["gene"]
                new.features.append(feat)
            else:
                continue
        # write the new features to a GFF
        GFF.write([new], out_handle)
        return
    out_handle.close()


if __name__ == "__main__":
    #GFF_bcftools_format(sys.stdin, sys.stdout)
    GFF_bcftools_format(sys.argv[1], sys.argv[2])
