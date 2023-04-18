import os
import sys

def find_best_template(res_file, database):
    """Returns the mapping results from the reference mapping"""
    if os.path.exists(res_file):
        template_score = 0
        reference_header_text = ""
        with open(res_file, 'r') as f:
            data = f.read().split("\n")
        data = data[:-1] #Last line is empty
        for item in data:
            item = item.split("\t")
            if item[0][0] != "#":
                if float(item[1]) > template_score:
                    template_score = float(item[1])
                    reference_header_text = item[0]
        template_number = findTemplateNumber(reference_header_text, database)
        return template_number, template_score, reference_header_text

def findTemplateNumber(name, database):
    if os.path.exists(database + ".name"):
        with open(database + ".name") as f:
            t = 1
            for line in f:
                if line.rstrip() == name:
                    return t
                else:
                    t += 1