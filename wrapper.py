#!/usr/local/anaconda/bin/python2.7
################################################################################
##########                    CGE Service Wrapper                    ###########
################################################################################
# This script is part of the CGE Service structure
# --> The input/arguments are validated
# --> Create service folder and subfolders
# --> The uploaded files are copied to the 'Upload' directory
# --> Log service submission in SQL
# --> Setup and execution of service specific programs
# --> The success of the service is validated
# --> A HTML output page is created if in webmode
# --> Downloadable files are copied/moved to the 'Download' Directory
# --> The SQL service entry is updated
# --> The files are zipped for long-term storage
import sys, os, time, random  # , ete3

# INCLUDING THE CGE MODULES (No need to change this part)
sys.path.append("/home/data1/services/CGEpipeline/CGEmodules")
from assembly_module import (PrepareAssembly, MakeAssembly,
                             printContigsInfo, AddAssemblyToProgList)
from functions_module_2 import (printDebug, copyFile, program, createServiceDirs,
                                getArguments, paths, makeFileList, fileUnzipper,
                                printInputFilesHtml, fileZipper, PrepareDownload,
                                UpdatePaths, proglist, CheckFileType, printOut,
                                moveFile, setDebug, add2DatabaseLog, dlButton,
                                GraceFullExit, FixPermissions, tsv2html)


################################################################################
##########                         FUNCTIONS                         ###########
################################################################################
def getIncPos(datafilesDir, fsaInput):
    datafiles = os.listdir(datafilesDir);
    skipLastChar = len("_alignment.fsa");

    printOut("<table style=\"width:100%\">\n");
    printOut("<tr bgcolor=#B22222><th>Isolate</th><th>Valid positions</th><th>Pct. of reference</th></tr>\n");

    for fsa in datafiles:
        if (fsa.endswith(".fsa")):
            name = fsa[:-skipLastChar];
            fileP = open(datafilesDir + fsa, "r");
            header = fileP.readline();
            inc = 0;
            seqlen = 0;
            for line in fileP:
                if (line.startswith(">")):
                    break;
                else:
                    seq = line.rstrip();
                    seqlen += len(seq);
                    inc += seq.count("A");
                    inc += seq.count("C");
                    inc += seq.count("G");
                    inc += seq.count("T");
                    if (fsaInput):
                        inc += seq.count("a");
                        inc += seq.count("c");
                        inc += seq.count("g");
                        inc += seq.count("t");
            fileP.close();

            printOut("<tr><td>%s</td><td>%d</td><td>%.2f</td><tr>\n" % (name, inc, 100.0 * inc / seqlen));
    printOut("</table><br>\n");

    return 0;


def makeIncOut(logFile, datafilesDir, fsaInput):
    inc = 0;
    tot = 1;
    fileP = open(logFile, "r");
    for line in fileP:
        line = line.rstrip();
        if (line.endswith("bases included in distance matrix.")):
            info = line.split(" ");
            inc = int(info[1]);
            tot = int(info[3]);
            break;
    fileP.close();

    printOut(
        "<b>Percentage of reference covered by all isolates: %.2f (%d / %d)</b><br>\n" % (100.0 * inc / tot, inc, tot));
    printOut(
        "Below is the single isolate stats on covered and trusted positions with respect to the reference.<br><br>\n");

    return getIncPos(datafilesDir, fsaInput);


def imitateUPGMA(filename):
    try:
        infile = open(filename, "r");
    except:
        sys.exit("No such file:\t%s\n", filename);

    newick = infile.readline();
    parts = newick.split(',');

    i = 1;
    num = len(parts);
    while (i < num):
        if ('(' not in parts[i] and ')' not in parts[i]):
            parts[0] = '(' + parts[0];
            parts[i] += '):0.00';
        i += 1;
    infile.close();
    newick = ','.join(parts);

    try:
        infile = open(filename, "w");
    except:
        sys.exit("Cannot open file:\t%s\n", filename);

    infile.write(newick);
    infile.close();

    return 0;


def printPhyloCanvas(paths, service, version, runID):
    # Make newick visible to the web
    newick = paths['downloads'] + 'results.nwck';
    tmpNewick = "/srv/www/htdocs/services/" + service + "-" + version + "/tmp/" + runID + ".nwck";
    os.system("sed \'s/\"//g\' %s > %s" % (newick, tmpNewick));
    FixPermissions(tmpNewick);
    
    # Make phylocanvas
    tmpCanvas = "/srv/www/htdocs/services/" + service + "-" + version + "/tmp/" + runID + ".php";
    outfile = open(tmpCanvas, "w");
    outfile.write("""<!DOCTYPE html>\n""");
    outfile.write("""<html>\n""");
    outfile.write("""<script type="application/javascript"\n""");
    outfile.write("""src="https://cge.food.dtu.dk/services/Evergreen-2.0/etc/phylocanvas-quickstart/phylocanvas-quickstart.js"></script>\n""");
    outfile.write("""<head>\n""");
    outfile.write("""\t<!--styling, obs. that .phylocanvas-history is hidden, could be in service css file-->\n""");
    outfile.write("""\t<style>\n""");
    outfile.write("""\tbody {\n""");
    outfile.write("""\t\tmargin: 0.625em auto;\n""");
    outfile.write("""\t\tmax-width: 80em;\n""");
    outfile.write("""\t}\n""");
    outfile.write("""\t#phylocanvas {\n""");
    outfile.write("""\t\twidth: 100%;\n""");
    outfile.write("""\t\theight: 500px;\n""");
    outfile.write("""\t}\n""");
    outfile.write("""\t.phylocanvas-history {\n""");
    outfile.write("""\t\tdisplay: none;\n""");
    outfile.write("""\t}\n""");
    outfile.write("""\t.results-content {\n""");
    outfile.write("""\t\twidth: 100%;\n""");
    outfile.write("""\t\theight: 500px;\n""");
    outfile.write("""\t\tbackground-color: #f0f0f0;\n""");
    outfile.write("""\t}\n""");
    outfile.write("""\t</style>\n""");
    outfile.write("""\t<body>\n""");
    # outfile.write("""\t\t<div class="results-content">\n""");
    # outfile.write("""\t\t\t<p></p>\n"""); # Text can go here
    # outfile.write("""\t\t</div>\n""");
    outfile.write("""\t\t<!--this div is used for the tree-->\n""");
    outfile.write("""\t\t<div id="phylocanvas"></div>\n""");
    # outfile.write("""\t\t<div class="results-content">\n""");
    # outfile.write("""\t\t\t<p></p>\n"""); # Some text below the tree
    # outfile.write("""\t\t</div>\n""");
    outfile.write("""\t\t<!--script for canvas object-->\n""");
    outfile.write("""\t\t<script type="application/javascript">\n""");
    outfile.write("""\t\tvar newick = <?php\n""");
    outfile.write("\t\t\t$string = file_get_contents(\"%s\");\n" % (runID + ".nwck"));
    outfile.write("""\t\t\t$string = str_replace("\n", "", $string);\n""");
    outfile.write("""\t\t\techo '"'.$string.'"';\n""");
    outfile.write("""\t\t?>;\n""");
    outfile.write("""\t\t//console.log(newick);\n""");
    outfile.write("""\t\tfunction defaultLoad() {\n""");
    outfile.write("""\t\t\tvar tree = Phylocanvas.createTree('phylocanvas');\n""");
    outfile.write("""\t\t\ttree.load(newick);\n""");
    outfile.write("""\t\t\ttree.setTreeType('rectangular');\n""");
    outfile.write("""\t\t\ttree.setNodeSize(4);\n""");
    outfile.write("""\t\t\ttree.setTextSize(24);\n""");
    outfile.write("""\t\t\ttree.lineWidth = 4;\n""");
    outfile.write("""\t\t\ttree.draw();\n""");
    outfile.write("""\t\t\treturn tree;\n""");
    outfile.write("""\t\t}\n""");
    outfile.write("""\t\tconst tree = defaultLoad();\n""");
    outfile.write("""\t</script>\n""");
    outfile.write("""</body>\n""");
    outfile.write("""\n""");
    outfile.close();
    FixPermissions(tmpCanvas);
    
    # Make inner html
    relNwck = "/services/" + service + "-" + version + "/tmp/" + runID + ".php";
    printOut("<iframe src=\"%s\" width=\"1024\" height=\"512\"><center><button type=\"button\" onclick=\"window.open('%s', '_self')\">PhyloCanvas</button></center></iframe>" % (relNwck, relNwck));
    
    return 0;


################################################################################
##########                           MAIN                            ###########
################################################################################
# SET GLOBAL VARIABLES
setDebug(False)
service, version = "MINTyper", "1.0"

# PARSE ARGUMENTS
# Add service specific arguments using the following format:
# (OPTION,   VARIABLE,  DEFAULT,  HELP)
# args = getArguments([
#   ('--uploadpath',  'uploadPath',  None, 'The folder containing uploads'),
#   ('-t',   'technology',  None, 'The sequencing platform of the input file')])
#
# Or by pasting the argument lines from the contigs file
args = getArguments('''
selectionbox   database            -d     VALUE  -d     ''
selectionbox   masking             -m     VALUE  -m     ''
text           prune               -p     VALUE  -p     ''
file           referenceFile       -ref   VALUE  -ref   ''
selectionbox   prunesignificance   -sig   VALUE  -sig   ''
text           cluster_length               -cluster_length     VALUE  -cluster_length     ''
''')

# VALIDATE REQUIRED ARGUMENTS
if args.database is None: GraceFullExit("Error: No database was chosen!\n")
if args.uploadPath is None or len(args.uploadPath) < 5 or not os.path.exists(args.uploadPath):
    GraceFullExit("Error: No valid upload path was provided! (%s)\n" % (args.uploadPath))
elif args.uploadPath[-1] != '/':
    args.uploadPath += '/'  # Add endslash

# SET RUN DIRECTORY (No need to change this part)
if args.reuseid:
    # /srv/www/secure-upload/isolates/5_16_4_2015_162_299_423438//0/
    runID = [x for x in args.uploadPath.split('/') if x != ''][-2]
else:
    runID = time.strftime('%w_%d_%m_%Y_%H%M%S_') + '%06d' % random.randint(0, 999999)
paths.serviceRoot = '{programRoot}IO/%s/' % (runID)
paths.isolateRoot = '{programRoot}'
paths.add_path('uploads', '{serviceRoot}uploads/')

# SET AND UPDATE PATHS (No need to change this part)
paths.add_path('webfolder', '/services/MINTyper-1.0/tmp/' + runID + '/');
paths.add_path('webdir', '/srv/www/htdocs' + paths['webfolder']);
UpdatePaths(service, version, '', '', args.webmode);

# CREATE SERVICE DIRECTORIES (No need to change this part)
createServiceDirs()
stat = paths.Create('uploads')

# LOG SERVICE SUBMISSION IN SQL (No need to change this part)
add2DatabaseLog(service + '-' + version, runID, args.usr, args.ip, args.database)

# MOVE UPLOADS FROM APP- TO ISOLATE UPLOAD DIRECTORY (No need to change this part)
if stat:
    # Move the uploaded files to the upload directory
    moveFile(args.uploadPath + '*', paths['uploads'])
    # Unzipping uploaded files if zipped
    # fileUnzipper(paths['uploads'])
    # GET INPUT FILES from input path
    inputFiles = makeFileList(paths['uploads'])
else:
    GraceFullExit("Error: Could not create upload directory!\n")

if len(inputFiles) > 100:
    GraceFullExit("Error: Attempt to analyze more than 100 files, please keep your uploads to 100 or less.\n")

# Get inputfiles splitted in Nanopore and non-nanopre sequences.
cmd = "%s/mintyper/fingerseq -i %s | grep -v \"^#\"> %s/fingerReport.tsv" % (
paths['scripts'], " ".join(inputFiles), paths['outputs']);
os.system(cmd);
infile = open(paths['outputs'] + "/fingerReport.tsv", "r");
Illumina_files = [];
Nano_files = [];
Asm_files = [];
PE = 0;
fsaInput = False;
for line in infile:
    line = line.rstrip();
    info = line.split("\t");
    if (info[1] == "Nanopore"):
        Nano_files.append(info[0]);
    elif(info[3] != "Na"):
        Illumina_files.append(info[0]);
    else:
        Illumina_files.append(info[0]);

maskingScheme = "";
if (args.masking == "DCM"):
    maskingScheme = "/home/data1/services/MINTyper/MINTyper-1.0/scripts/mintyper/dcmMethylations.txt";

# ADDING PROGRAMS
mintyper = proglist.AddProgram(program(
    name='MINTyper-1.1.0', path=paths['scripts'] + 'mintyper/bin/mintyper', timer=0,
    ptype='python3', toQueue=True, wait=False, workDir='',
    args=['--prune_distance', args.prune,
          '--output', paths['outputs']]))

# dummy.AppendArgs(['-option', 'content'])
# here
# Remove when prev options are set back
# mintyper.AppendArgs(["-i_path", paths['uploads']]);
# Remove False statements when prev options are set back
if (len(Nano_files) != 0):
    mintyper.AppendArgs(["--nanopore"] + Nano_files);
    #input_path = '{}/*'.format(paths['uploads']);
    #mintyper.AppendArgs(['--nanopore', input_path]);
if (len(Illumina_files) != 0):
    mintyper.AppendArgs(["--illumina"] + Illumina_files);
    #input_path = '{}/*'.format(paths['uploads']);
    #mintyper.AppendArgs(['--illumina', input_path]);

if (len(maskingScheme) != 0):
    mintyper.AppendArgs(['--masking_scheme', maskingScheme]);

if (int(args.cluster_length) > 0):
    mintyper.AppendArgs(['--cluster_length', args.cluster_length]);

if (args.prunesignificance == "insignificant"):
    mintyper.AppendArgs(['-insig_prune']);

# MISSING DB PATH
dbPath_host = "/home/data2/databases/internal/kma/kmerdb/";
reffilename = "";
if (args.referenceFile != None):
    # dbPath_host += "test_kma-1.0.1/bacteria/bacteria.ATG";
    moveFile(args.referenceFile, paths['uploads']);
    reffilename = paths['uploads'] + args.referenceFile.split('/')[-1];
    dbPath_host = "";
elif (args.database == "bacteria"):
    dbPath_host += "test_kma-1.0.1/bacteria/bacteria.ATG";

if (dbPath_host != ""):
    mintyper.AppendArgs(['--db', dbPath_host]);
else:
    mintyper.AppendArgs(['--reference', reffilename]);

# EXECUTION OF THE PROGRAMS
mintyper.Execute(True)
mintyper.WaitOn(interval=10)

# THE SUCCESS OF THE SERVICE IS VALIDATED
status = mintyper.GetStatus()
if status != 'Done': GraceFullExit("Error: Execution of the program failed!\n")

# CREATE THE HTML OUTPUT FILE (No need to change this part)
if args.webmode:
    if os.path.exists(paths['outputs']) and not os.path.exists(paths['outputs'] + service + ".out"):
        with open(paths['outputs'] + service + ".out", 'w') as f:
            f.write('<h1>%s-%s Server - Results</h1>' % (service, version))

# PRINT THE STANDARD OUTPUT OF THE PROGRAM TO THE HTML OUTPUT FILE
# mintyper.Printstdout()

# Print warning if only a single isolate was uploaded
if ((len(Illumina_files) + len(Nano_files)) == 1):
    printOut("<font color=\"#B22222\" size=\"6\"><b>WARNING</b><br>\n");
    printOut("Please note that at least two isolates are needed to provide a distance between them.</font><br><br>\n");

# Make tree visualization
# os.system("%s %s%s %s" %("/home/data1/tools/bin/Anaconda3-2.5.0/bin/python3", paths['scripts']+'bin/', "nwck2png.py", paths['serviceRoot']+'outputs/outtree.newick'));
'''
treepath = paths['serviceRoot']+'outputs/outtree.newick';
try:
  tree = Tree(treepath)
except NewickError as e:
  exiting("Couldn't open {0}".format(treepath))
graphpath = os.path.join(paths['serviceRoot']+'outputs/', "outtree.PNG")
style = TreeStyle()
style.show_branch_length = True
style.show_leaf_name = True
#style.branch_vertical_margin = 10
tree.render(graphpath, tree_style=style, dpi=300)
'''
# Print tree to output website
# printOut("<br>");
# printOut("<img src=\"%s\" width=\"950\" height=\"400\" border=\"1\">" %(paths['serviceRoot']+'outputs/cgeout/outtree.PNG'));
# printOut("<p><div class=\"bulk\" align=\"justify\">");

# Make Zip archive with vcf files
cmd = "mkdir %s/results_vcf/" % (paths['outputs']);
os.system(cmd);
cmd = "for vcf in %s/alignments/*.vcf.gz; do mv $vcf ${vcf%%.vcf.gz}.vcf.txt.gz; mv ${vcf%%.vcf.gz}.vcf.txt.gz %s; done" % (
paths['outputs'], paths['outputs'] + 'results_vcf/');
os.system(cmd);
cmd = "gunzip %s/*" % (paths['outputs'] + 'results_vcf');
os.system(cmd);
cmd = "zip -j -r -q %s/results_vcf %s/results_vcf/ && rm -rf %s/vcf/*" % (
paths['outputs'], paths['outputs'], paths['outputs']);
os.system(cmd);

# cmd = "zip -r %s %s && rm -rf %s" %(paths['downloads'] + 'results_vcf', paths['outputs'] + 'vcf/', paths['outputs'] + 'vcf/');
# os.system(cmd);

# Mv files to download
moveFile(paths['serviceRoot'] + 'outputs/mintyper.log', paths['downloads'] + 'results.log');
moveFile(paths['serviceRoot'] + 'outputs/cluster.dbscan', paths['downloads'] + 'cluster.dbscan');
moveFile(paths['serviceRoot'] + 'outputs/distmatrix.txt', paths['downloads'] + 'results.txt');
moveFile(paths['serviceRoot'] + 'outputs/tree.newick', paths['downloads'] + 'results.nwck');
moveFile(paths['serviceRoot'] + 'outputs/results_vcf.zip', paths['downloads'] + 'results_vcf.zip');
#moveFile(paths['serviceRoot'] + 'outputs/template_sequence*', paths['downloads'] + 'template_sequence.fasta');

# Make phylocanvas
if os.path.exists(paths['downloads'] + 'results.nwck'):
    printPhyloCanvas(paths, service, version, runID);

# Make TreeViewer Link, Deprecated
if False and sos.path.exists(paths['downloads'] + 'results.nwck'):
    paths.Create('webdir');
    copyFile(paths['downloads'] + 'results.nwck', paths['webdir'] + 'results.nwck');
    # TreeViewer cannot handle NJ trees, so we imitate UPGMA here
    # imitateUPGMA(paths['webdir'] + 'results.nwck');
    printOut("""<svg id='svg'
  xmlns='http://www.w3.org/2000/svg'
  xmlns:xlink='http://www.w3.org/1999/xlink'
  version='1.1'></svg><br>""")
    printOut("""
  <button onclick='svgdownload("png");'>Download as PNG</button>
  <button onclick='svgdownload("pdf");'>Download as PDF</button>""")
    printOut("""
  <script type='text/javascript'>
    Newick2Tree('%s/%s', [1024, 512], 360, 0, 0, ['8', '000000', 'Verdana, Geneva, sans-serif'], ['6', '000000', 'Verdana, Geneva, sans-serif'], false);
    d3.selectAll('circle.node').style('display', 'None');
    $(document).ready(function(){
      d3.selectAll('circle.node').style('display', 'None');
    });
  </script>""" % (paths['webfolder'], 'results.nwck'));
    treeviewerLink = '/services/TreeViewer-1.0/index.php?imagedir=%s&newick=%s' % (paths['webfolder'], 'results.nwck');
    printOut('Click <a href="%s">here</a> to modify the tree in the advanced TreeViewer.<br><br>\n' % (treeviewerLink));

# PREPARE THE DOWNLOADABLE FILE(S)
# _ = PrepareDownload(paths['serviceRoot']+'outputs/results.log')
# _ = PrepareDownload(paths['serviceRoot']+'outputs/results.txt')
# _ = PrepareDownload(paths['serviceRoot']+'outputs/results.nwck')
# _ = PrepareDownload(paths['serviceRoot']+'outputs/results.aln')

# PRINT STATS
makeIncOut(paths['downloads'] + 'results.log', paths['serviceRoot'] + 'outputs/alignments/', fsaInput);

# PRINT TAB FILE(S) AS RESULT TABLES TO THE HTML OUTPUT FILE
# tsv1=paths['downloads']+'results.res'

# PRINT THE PROGRAM PARAMETERS TO THE HTML OUTPUT FILE
# if len(inputFiles) > 0: printInputFilesHtml(inputFiles)

# PRINT FILE DOWNLOAD-BUTTON(S) TO THE HTML OUTPUT FILE
dlButton('Log', 'results.log')
dlButton('Distance matrix', 'results.txt')
dlButton('Phylogentic tree', 'results.nwck')
dlButton('Vcf files of mutations', 'results_vcf.zip')
#dlButton('Reference Sequence', 'template_sequence.fasta')
dlButton('Cluster.dbscan', 'cluster.dbscan')

# dlButton('Gene mapping results', 'results.res')

## PRINT THE EXTENDED OUTPUT
# aln = paths['downloads']+'results.aln'
# if(os.path.exists(aln)):
#   printOut('<br>\n')
#   printExtendOutput(aln)
#   printOut('<br>\n')

################################################################################
# LOG THE TIMERS (No need to change this part)
proglist.PrintTimers()

# FIX THE PERMISSIONS OF THE SERVICE ROOT
FixPermissions()

# INDICATE THAT THE WRAPPER HAS FINISHED (No need to change this part)
printDebug("Done")

# GZIP ALL FILES IN THE SERVICE DIRECTORY AND SUBDIRS FOR LONG-TERM STORAGE (No need to change this part)
fileZipper(paths['outputs'])
# fileZipper(paths['serviceRoot'])


