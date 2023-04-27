import os
import sys

def find_best_template_from_spa_file(spa_file, database):
    """Returns the mapping results from the reference mapping"""
    if os.path.exists(spa_file):
        template_score = 0
        template_number = None
        reference_header_text = ""
        with open(spa_file, 'r') as f:
            data = f.read().split("\n")
        data = data[:-1] #Last line is empty
        for item in data:
            item = item.split("\t")
            if item[0][0] != "#":
                if float(item[2]) > template_score:
                    template_score = float(item[2])
                    template_number = int(item[1])
                    reference_header_text = item[0]
        return template_number, template_score, reference_header_text
