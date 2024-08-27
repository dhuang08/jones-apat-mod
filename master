#module imports
import json
import urllib.request
from pip._vendor import requests

#pip3 install biopython   <- required
from Bio import SeqIO
#Function for transforiming decimal scores to colour ranges (direct translation)
def cgrad(numin):
    hue = float(numin)*0.7
    hue = 0.7-hue    
    sixth = int(6.0 * hue)
    rising  = (hue - (sixth / 6.0)) * 6.0
    falling = 1.0 - rising
    InvSat  = 1.0 - 1.0
    red =0
    green = 0
    blue =0
    if((sixth == 0) or (sixth == 6)):
    
        red   = 1.0
        green = rising
        blue  = 0.0
    
    elif sixth == 1:
    
        red   = falling
        green = 1.0
        blue  = 0.0
    
    elif(sixth == 2):
    
        red   = 0.0
        green = 1.0
        blue  = rising
    
    elif(sixth == 3):
    
        red   = 0.0
        green = falling
        blue  = 1.0
    
    elif(sixth == 4):
    
        red   = rising
        green = 0.0
        blue  = 1.0
    
    elif(sixth == 5):
    
        red   = 1.0
        green = 0.0
        blue  = falling
    


    red   += ((1.0-(red))   * InvSat)
    green += ((1.0-(green)) * InvSat)
    blue  += ((1.0-(blue))  * InvSat)

    red   *= 255
    green *= 255
    blue  *= 255

    red   = "%02lx"% int(red)
    green = "%02lx"% int(green)
    blue  = "%02lx"% int(blue)
    
    return("#%s%s%s"% (red,green,blue))

#file opens,input.txt is for input in gff or csv, htmltest1 is the output html, out.csv is output csv
inf = open("input1.txt", "r")
file = open("htmltest1.html", "wt")
csfile = open("out.csv", "wt")

#global variables, inlist is for reading input1.txt (each line to list)
#pred list has all the prediction values, ex 0.597
inlist = []
predlist= []
for i in inf:
    inlist.append(i)
#seq is protein sequence
seq = ""
link = ""
tempc = 0
seqid = ""

# **only seq and predlist are required for proper output 

def ingff():
    global seq
    global link
    global predlist
    global tempc

    for i in range(6, inlist.index("##end-Protein\n")):
        temps = inlist[i]
        seq += temps[2:len(temps)-1]
    seq = seq.strip()
    predlist= ["0.000"]*len(seq)
    tempc = len(seq)%60
    prelink = inlist[0].split()
    link += "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?"+ prelink[0]
    tlines = int((len(seq) - tempc)/60 +1)
    tint = 0
    for i in range(9+tlines, len(inlist)):
        if(inlist[i].find("phos-unsp")!=-1):
            templ = inlist[i].split()
            tint = int(templ[4])
            predlist[tint-1] =templ[5]

def incsv():
    global seq
    #global link
    global predlist
    global tempc

    #still needs more work : exact specification of input file required
    #right now, sees as no headers, first column seq letter, next column value.
    ci = 0
    for i in inlist:
        ts = i.split(",")
        seq += ts[0]
        predlist.append("0.000")
        predlist[ci] = ts[1]

def injsn():
    global predlist
    global seq
    global tempc
    global seqid
    ffile = open("infastemp.fasta", "+r")
    
    with open('jsin.json', encoding='utf-8') as fh:
        tdata = json.load(fh)
    sl = ("https://rest.uniprot.org/uniprotkb/%s.fasta"% str(tdata['options']['accession']).strip("['").strip("']"))
    
    try:
        with urllib.request.urlopen(sl) as f:
            
            ffile.write(f.read().decode('utf-8'))
            ws = SeqIO.parse("infastemp.fasta", 'fasta')
            tnext = next(ws)
            seqid = tnext.description
            seq = tnext.seq
            predlist = ["0.000"]*len(seq)
    except urllib.error.URLError as e:
        print("Accessing uniprot.org Error: ", e.reason)
    
    for i in range(0, len(tdata['data'])):
        predlist[i] = (tdata['data'][str(i+1)])
        #print(tdata['data'][str(i+1)])

    ffile.close()


def html():
    global tempc
    text1 = """
<html>
<head>
    <title>Automated Protein Annotation Tool</title>
    <style type="text/css">
    <!--
    .header {font: bold 32pt Arial,Helvetica,sans-serif;
             text-align: center;
             color: orange;
             padding: 2px 0px 1em 0px;
            }
    .highlight {
                background: orange;
               }
    .warning {
              background: red;
              font: 18pt Arial,Helvetica,sans-serif;
              text-align: center;
             }
    .results {
              background: #ffffff;
              border: thin solid black;
              color: #000000;
              margin: 20px;
              padding: 0px 0px 10px 0px;
             }
    h1 {font: 18pt Arial,Helvetica,sans-serif;
        text-align: left;
        color: white; 
        background: #336699;
        padding: 4px;
        margin: 0px;
       }
    h2 { padding: 2px 4px 4px 4px;
         margin: 0px;
       }
    h3 { padding: 2px 4px 4px 10px;
         margin: 0px;
       }
    h4 { padding: 2px 4px 4px 10px;
         margin: 0px;
       }
    p  { padding: 0px 4px 4px 24px;
         margin: 0px;
       }
    table { margin: 2px 10px 2px 10px;
       }
    -->
    </style>
</head>
<body bgcolor="white">
<p class="header">Automated Protein Annotation Tool</p>
<div class='results'>
<h1> INPUT DETAILS</h1>
<p><b>Number of Residues:</b> %d </p>
<p><b>Sequence ID:</b> %s</p>
<p><b>Sequence:</b></p>
<pre><p> """ % (len(seq), seqid)
    file.write(text1)

    for i in range(0, len(seq)-60, 60):
        file.write("\t")
        for j in range(i, i+60):
            file.write(seq[j])
    
        file.write("<br />")
    tempc = len(seq)%60
    if(tempc!=0):
        file.write("\t")
        for i in range(len(seq)-tempc, len(seq)):
            file.write(seq[i])
        file.write("<br />")

    file.write("""</p></pre>
 <p><b>Parameters:</b></p>
</div>
<hr /> <div class='results'>
<h1>Index of tools run:</h1><br /><i>Click on the name of the tool to avoid scrolling down.</i><br />
<a href="#NetPhos">
    <h3>NetPhos</h3>
</a>
</div>
<hr />
""")

    file.write("""
<div class='results'>
<a name= NetPhos ><h1>Analysis results from program: NetPhos </h1></a>
<h2>Running NetPhos Version 3.1</h2>
<h3>Function : Protein Phosphorylation Sites Prediction - Unspecified/Overall Score</h3>
<h3><a href='https://services.healthtech.dtu.dk/services/NetPhos-3.1/'>NetPhos Web Server</a></h3>    
<h3><a href='%s'>Actual prediction(native, unparsed form)- available only for a limited time</a></h3>
""" % link)

    for i in range(0, len(seq)-15, 15): 
        file.write("""<p></p><table border="1">
    <tbody>
    <tr><td colspan='2'>Residue</td>""")
        for j in range(i, i+15):
        
            if(float(predlist[j])>0.5):
                file.write("""<td width='30' bgcolor='red'>%s</td>"""% seq[j])
            else:
                file.write("""<td width='30'>%s</td>"""% seq[j])
        file.write("""</tr>
    <tr><td colspan='2'>Number</td>
    """)
        for j in range(i, i+15):
            file.write("""<td width='30'>%d</td>"""% int(j+1))
        file.write("""</tr
<tr><td rowspan='1'>NetPhos</td>
<td>P-Score</td>
""")
        for j in range(i, i+15):
            file.write("""<td width='30' bgcolor='%s'>%0.5s</td>"""% (cgrad(predlist[j]), predlist[j]))
        file.write("""</tr>
    </table>""")

    if(len(seq)%15!=0):
        file.write("""<p></p><table border="1">
    <tbody>
    <tr><td colspan='2'>Residue</td>""")
        for j in range(len(seq)-len(seq)%15, len(seq)):
            if(float(predlist[j])>0.5):
                file.write("""<td width='30' bgcolor='red'>%s</td>"""% seq[j])
            else:
                file.write("""<td width='30'>%s</td>"""% seq[j])
        file.write("""</tr>
    <tr><td colspan='2'>Number</td>
    """)
        for j in range(len(seq)-len(seq)%15, len(seq)):
        
            file.write("""<td width='30'>%d</td>"""% int(j+1))
        file.write("""</tr
<tr><td rowspan='1'>NetPhos</td>
<td>P-Score</td>
""")
        for j in range(len(seq)-len(seq)%15, len(seq)):
            file.write("""<td width='30' bgcolor='%s'>%0.5s</td>"""% (cgrad(predlist[j]), predlist[j]))
        file.write("""</tr>
    </table>""")
    file.write("</tbody>")
    file.write(""" 
</div>
<p><br /><br /></p><hr /><p><br /><br /></p>
              
</body>
</html>
""")

def csv():
    tcount = 1
    csfile.write("track,t_name,t_type,t_position,t_colour,t_opacity,t_text_colour,t_help,entry,text,value,sequence,hover,link,start,end,position,colour,opacity,text_colour")
    for i in predlist:
    
        csfile.write("\n2,TEST2,histogram,-1,,0.5,,test track histogram,%d,,%s,%s, %s (%s),,,,%d,%s,0.9," % (tcount, i, seq[tcount-1], seq[tcount-1], i, tcount, cgrad(i)))
        tcount+=1

uin = input("Input.txt or jsin.json? [t/j]\n")
match uin:
    case "t" |"T":
        txtin = input("GFF or CSV? [g/c]\n")
        match txtin:
            case "G" | "g":
                ingff()
                print("Run Successful")
            case "C" | "c":
                incsv()
                print("Run Successful")
            case _:
                print("Invalid input, run again")
    case "j" | "J":
        injsn()
        print("Run Successful")
    case _:
        print("Invalid input, run again")
html()
csv()





csfile.close()
inf.close()
file.close()
