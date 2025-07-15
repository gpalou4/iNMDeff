"Modified functions from bcbio-nextgen"

import os
import sys

def postprocess_somatic(in_file):
    """Post-process somatic calls to provide standard output.
    - Converts SGT and NT into standard VCF GT fields
    - Replace generic TUMOR NORMAL names in VCF with sample names.
    """
    out_file = in_file.replace(".vcf", "_GT_fixed.vcf")
    with open(in_file, "r") as in_handle:
        with open(out_file, "w") as out_handle:
            added_gt = False
            normal_index, tumor_index = (None, None)
            for line in in_handle:
                if line.startswith("##FORMAT") and not added_gt:
                    added_gt = True
                    out_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    out_handle.write(line)
                elif line.startswith("#CHROM"):
                    assert added_gt
                    parts = line.strip().split("\t")
                    normal_index = parts.index("NORMAL")
                    tumor_index = parts.index("TUMOR")
                    #line = line.replace("NORMAL", paired.normal_name).replace("TUMOR", paired.tumor_name)
                    out_handle.write(line)
                elif line.startswith("#"):
                    out_handle.write(line)
                else:
                    parts = line.rstrip().split("\t")
                    tumor_gt, normal_gt = tumor_normal_genotypes(parts[3], parts[4].split(","),
                                                                    parts[7].split(";"), in_file, parts[:2])
                    parts[8] = "GT:%s" % parts[8]
                    parts[normal_index] = "%s:%s" % (normal_gt, parts[normal_index])
                    parts[tumor_index] = "%s:%s" % (tumor_gt, parts[tumor_index])
                    out_handle.write("\t".join(parts) + "\n")
        return out_file

def tumor_normal_genotypes(ref, alt, info, fname, coords):
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
    """
    known_names = set(["het", "hom", "ref", "conflict"])
    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "confict"]):
            return "0/0"
        else:
            # Non-standard representations, het is our best imperfect representation
            # print(fname, coords, ref, alt, info, val)
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref] + alt)}
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt
    nt_val = [x.split("=")[-1] for x in info if x.startswith("NT=")][0]
    normal_gt = name_to_gt(nt_val)
    sgt_val = [x.split("=")[-1] for x in info if x.startswith("SGT=")]
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val[0].split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return tumor_gt, normal_gt

in_file = sys.argv[1]
out_file = postprocess_somatic(in_file)