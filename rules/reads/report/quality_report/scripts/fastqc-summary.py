#! python

from argparse import ArgumentParser
import os

fastqc_suffix = '_fastqc.html'

parser = ArgumentParser(
    description="Takes all fastqc reports in the specified directory and generated one page summarized report")
parser.add_argument('summary_file', help="Output html file with summarized fastq files")
parser.add_argument('reports', nargs="*", help="Fastqc files (files with %s suffix)" % fastqc_suffix)
args = parser.parse_args()

report_files, out_file = sorted(args.reports), args.summary_file

from xml.dom import minidom

report_names = [os.path.basename(rf)[:-len(fastqc_suffix)] for rf in report_files]
reports = [minidom.parse(open(rf)) for rf in report_files]


def retrieve_headers(report):
    return report.getElementsByTagName('h2')


def retrieve_icons(report):
    h2s = retrieve_headers(report)[1:]
    return map(lambda x: x.getElementsByTagName('img')[0], h2s)


def retrieve_labels(report):
    h2s = retrieve_headers(report)[1:]
    return map(lambda x: x.lastChild.toxml(), h2s)


def retrieve_summary(report):
    stat_table = report.getElementsByTagName('table')[0].childNodes[1]
    sequences_row, lens_row = 3, 5
    return map(lambda row_id: stat_table.childNodes[row_id].childNodes[1].firstChild.nodeValue,
               [sequences_row, lens_row])


icons = map(retrieve_icons, reports)
labels = retrieve_labels(reports[0])

doc = minidom.Document()
html, head, body, table = map(lambda x: doc.createElement(x), ['html', 'head', 'body', 'table'])
doc.appendChild(html)
html.appendChild(head)
html.appendChild(body)
body.appendChild(table)


def prepare_text_cell(text):
    td = doc.createElement('td')
    td.attributes['style'] = 'padding-left: 20px; padding-right: 20px'
    td.appendChild(doc.createTextNode(text))
    return td


def generate_icon_row(report_name, report_file, report_icons, report):
    tr = doc.createElement('tr')
    sequences, lens = retrieve_summary(report)
    for text in [report_name, '{0:,}'.format(int(sequences)), lens]:
        tr.appendChild(prepare_text_cell(text))
    for icon, pos in zip(report_icons, range(len(list(report_icons)))):
        td, a = doc.createElement('td'), doc.createElement('a')
        a.attributes['href'] = '%s#M%s' % (report_file.split('/')[-1], pos)
        a.attributes['title'] = labels[pos]
        a.appendChild(icon)
        td.appendChild(a)
        tr.appendChild(td)
    return tr


rows = []
for (report_name, report_file, report_icons, report) in zip(report_names, report_files, icons, reports):
    rows.append(generate_icon_row(report_name, report_file, report_icons, report))

for row in rows:
    table.appendChild(row)

doc.writexml(open(out_file, 'w'))
