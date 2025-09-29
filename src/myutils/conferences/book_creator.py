# Code used to create the book of abstracts and the program.
# Initially created by Daniel Sucerquia (dsucerg@gmail.com) for the
# Simplaix 2025 workshop.

import pandas as pd
from myutils.miscellaneous import output_terminal
import unicodedata


class AbstractBook:
    def __init__(self):
        self.order_talks = []
        self.order_posters = {}

        # days of the week from short to long
        self.short2long = {"Mon": "Monday",
                           "Tue": "Tueday",
                           "Wed": "Wednesday",
                           "Thu": "Thursday",
                           "Fri": "Friday",
                           "Sat": "Saturday",
                           "Sun": "Sunday"}
        
        # words to be avoided to be modified as title format
        self.exceptions = {"a", "an", "and", "as", "at", "but", "by", "for",
                           "in", "nor", "of", "on", "or", "so", "the", "to",
                           "up", "yet", "with", "through", "into", "is"}
        
        # add acronyms to avoid to transform them in title format
        self.acronyms = {}

        self.setter = {'selector': '', # filters Nan contribution of this kind
                       "contri-title" : "", # title of the contribution
                       "contri-abstract" : ""} #abstract of the contribution

    # ==== General ============================================================
    def section_separation(self, title):
        """
        Adds a new page separating the sections.

        Parameters
        ==========
        title: str
            title of the section.
        
        Return
        ======
        (str) A new page separating sections with a big title in an isolated
        page.
        """
        text = f"""

\\newpage
\\thispagestyle{{empty}}
\\vspace*{{\\fill}}
\\begin{{center}}
    \\Huge \\color{{hitsblue}}{{ \\textbf{{ {title} }} }}
\\end{{center}}
\\vspace*{{\\fill}}
\\newpage
"""
        return text


    def simplaix_main(self, subfiles, output):
        text = """
\\documentclass{article}
\\usepackage{xcolor}
\\usepackage{setspace}
\\usepackage{longtable}
\\usepackage{array}
\\usepackage{graphicx} % Required for inserting images
\\usepackage{pdfpages}
\\usepackage[utf8]{inputenc}
\\usepackage{imakeidx}
\\makeindex[name=auth, title={\centering\color{hitsblue}Author Index}, program=makeindex, options=-s authorstyle.ist]

\\definecolor{hitsblue}{HTML}{1f4a87}
\\newcommand{\\authorentry}[2]{\\index[auth]{#1@#1, #2}}

\\title{abstract book simplaix}
\\author{dsucerg }
\\date{April 2025}

\\def\\Ntag{%
  \\leavevmode
  \\hbox{%
    $\\rm N^{\\mkern0.8mu\\underline{\\mkern-0.8mu o\\mkern-0.8mu}\\mkern0.8mu}$%
  } \\hspace{0mm}
}

\\begin{document}

\\includepdf[pages=-]{cover.pdf}
"""

        for subfile in subfiles.split(','):
            text += f"\\input{{ {subfile} }}\n"
        
        text += "\n\\printindex[auth]\n\\end{document}\n"

        with open(output, 'w') as outfile:
            outfile.write(text)
        
        return text

    # ==== Program ================================================================
    def add_date(self, day):
        """
        Add the date of the day before before starting to show the events.

        Parameters
        ==========
        day: str
            date of the day
        
        Return
        ======
        (str) text in LaTeX format where a new page is created starting with the
        date as a title.
        """
        text = f"""

\\newpage

{{\\setlength{{\\parindent}}{{0pt}}

\\textbf{{\\Large \\color{{hitsblue}}{{ {day} }} }}
\\vspace{{1 cm}}

    """
        return text

    def session_title(self, event):
        """
        Add the title of the session including the chair before showing the
        events of the session.

        Parameters
        ==========
        event: dict
            Title of the session in event-dictionary format with all the
            information in the "title" key.
        
        Return
        ======
        (str) text in LaTeX format where the centered title of the session is
        added.
        """
        session, chair = event["title"].split(" &#8211; Chair: ")
        text = f"""
\\begin{{center}}
\\textbf{{\\color{{hitsblue}}{{ {session} }} }} \\\\{{\\bf Chair:}} {chair}
\\end{{center}}"""
        return text

    def add_item(self, event, timespace=1, titlespace=10.5):
        time = event["startTime"][17:22]
        title = event["title"].replace("&#8211;", "--").replace("\\/", "/")
        # Those with authors, add page reference
        if ']' in title and '--' in title:
            author_name = title.split('--')[0]
            author_name = author_name.split(']')[1].lower().replace(" ", "")
            author_name = self.remove_tildes(author_name)
            label = f" -- (page \\pageref{{talks:{author_name}}})"
            self.order_talks.append(author_name)
        else:
            label = ""
        text = f"""
\\begin{{tabular}}{{p{{{timespace}cm}}p{{{titlespace}cm}}}}
    {time} & {title}{label} \\\\
\\end{{tabular}}
    """
        return text

    def simplaix_program(self, link, output):
        """
        Copy the schedule from the webpage and create the program.

        Parameters
        ==========
        link: str
            link to the webpage to the Simplaix workshop.

        Return
        ======
        (str) Program in LaTeX format.
        """
        # Import events
        all_events = output_terminal(f'curl {link} ' +\
                                     '| grep \'"items"\' | ' +
                                     'sed "s/<\/script>//g" ',
                                    print_output=False)
        all_events = eval(all_events)

        # separate events by the day it happens
        events_per_day = {}
        for event in all_events["items"]:
            if event["startTime"][:16] not in list(events_per_day.keys()):
                events_per_day[event["startTime"][:16]] = []
            events_per_day[event["startTime"][:16]].append(event)

        # creates text
        text = self.section_separation('Program') # Page separating

        for day in events_per_day:
            text += self.add_date(day.replace(day[:3],
                                              self.short2long[day[:3]]))

            session = False
            for event in events_per_day[day]:
                if "Session" ==  event["title"][:7]:
                    if "&#8211; Chair:" in  event["title"]:
                        text += self.session_title(event)
                        session=True
                    else:
                        title_and_kind = event["title"].split("\\/")[1]
                        kind = title_and_kind.split(":")[0]
                        author = title_and_kind.replace(kind + ": ", "")
                        if event['excerpt'] != '':
                            title = event['excerpt'].replace("&#8220;",
                                                             "``")
                            title = title.replace("&#8221;", "''")
                        event["title"] = f" [{kind}] {author} -- {title}"
                        text += self.add_item(event)

                        continue
                else:
                    if session:
                        text += "\n\\vspace{0.25 cm}\n"
                        session = False
                    text += "\n\\vspace{0.25 cm}\n"
                    text += self.add_item(event)
            text += "\n}\n"

        with open(output, 'w') as outfile:
            outfile.write(text)

        return text


    # ==== Abstracts ==============================================================
    def costume_title(self, text):
        """
        Transform to title format.

        Parameters
        ==========
        text: str
            title as provided by the authors.
        
        Return
        ======
        (str) modified text in title style.

        Note
        ----
        Inside of the script, there are exceptions of articles (a, the, and ..) and
        words that should not be transformed such as self.acronyms. Access and modify
        according to your preferences. 
        """
        
        words = text.split()
        titled = [words[0] if words[0] in self.acronyms else
                  words[0].capitalize()]  # Always capitalize the first word
        for word in words[1:]:
            word = word if ((word in self.exceptions) or (word in self.acronyms)) \
                             else word.capitalize()
            # capitalize also after -
            if '-' in word and word not in self.acronyms:
                start = [word[0]]
                word = ''.join(start + [word[i].upper() if word[i - 1] == '-'
                                        else word[i]
                                        for i in range(1, len(word))])
            titled.append(word)
        return ' '.join(titled)

    def title_contri(self, title, label=''):
        """
        Title of the contribution which could be a talk, a poster or whatever
        with title, author and abstract.
        
        Parameters
        ==========
        title: str
            title of the contribution.

        Return
        ======
        (str) text in latex format for the title of the contribution.
        """
        text =f"""

\\begin{{center}}
    \\textbf{{\\Large \\color{{hitsblue}}{{ {self.costume_title(title)} }}  }}
    {label}
\\end{{center}}
"""
        return text


    def author_line(self, authors_info, affiliation):
        """
        Add author name, and affiliation.

        Parameters
        ==========
        authors_info: str
            name and e-mail.
        affiliation: str
            institution, city, country.
        
        Return
        ======
        (str) text in LaTeX formatwith author information adding some space
        before and after.
        """
        text =f"""

\\begin{{center}}
    \\vspace{{0.5cm}}
    \\textbf{{\\bf{{ {authors_info} }} }}\\\\
    {affiliation}
    \\vspace{{0.5cm}}
\\end{{center}}

"""
        return text

    def contribution(self, row, kind, stand):
        """
        Create the whole page of the contribution.

        Parameters
        ==========
        row: DataFrame
            Information about the contribution.
        
        Return
        ======
        (str) text in LaTeX format of the contribution.
        """
        name_label = f"{kind}:" + row["firstName"].replace(" ", "") + \
                    row["lastName"].replace(" ", "")
        name_label = self.remove_tildes(name_label).lower()

        title = self.costume_title(row[self.setter["contri-title"]])
        label = f'\\authorentry{{{ row["lastName"].title() }}}' + \
                f'{{{ row["firstName"].title() }}}' + \
                f'\\label{{{name_label}}}'
        text = self.title_contri(title, label=label)

        name = row["firstName"].title()
        lastname = row["lastName"].title()
        email = row["email"]
        complete_name = name + " " + lastname + " (" + email + ')'
        accepted_titles = {'Dr.': 'Dr. ',
                           'Dr': 'Dr. ',
                           'Prof. Dr.': 'Prof. Dr. ',
                           'Professor': 'Prof. ',
                           'Prof.': 'Prof. '}
        
        if row['title'] in accepted_titles.keys():
            complete_name = accepted_titles[row['title']] + complete_name
        
        if self.toc:
            self.order_posters[name_label] = [stand,
                                              title,
                                              name + " " + lastname]

        affiliation = row["institute"] + ', ' + row["city"] + ', ' + \
                      row['country']
        text += self.author_line(complete_name, affiliation)

        text += row[self.setter["contri-abstract"]] + '\n\n'
        if kind == 'Posters':
            text += '\\begin{flushright}\n' + \
                    f'    Stand \\Ntag {stand}\n' + \
                    '\\end{flushright}\n'
        text += '\\newpage\n\n'

        return text

    def remove_tildes(self, text):
        normalized = unicodedata.normalize('NFD', text)
        return ''.join(c for c in normalized if unicodedata.category(c) != 'Mn')

    def reorder_talks(self, subset):
        ordered = []
        for target_name in self.order_talks:
            found = False
            for _, row in subset.iterrows(): 
                name_label = row["firstName"].replace(" ", "") + \
                             row["lastName"].replace(" ", "")
                name_label = self.remove_tildes(name_label)
                
                if name_label.lower() == target_name:
                    found = True
                    ordered.append(row)
                    break

            if not found:
                print(f"{target_name} was not found in the abstracts")
        return ordered

    # add2executable
    def abstracts(self, section, excel_file, output, reorder=False,
                  toc=False):
        """
        Write the whole abstract book in a file, usually a .tex because it is
        LaTeX formated.

        Parameters
        ==========
        excel_file: str, path
            excel file with the information of the abstracts, titles and others.
        output: str, path
            output file.
        
        Parameters
        ==========
        (str) text added to the file.
        """
        self.toc = toc
        df = pd.read_excel(excel_file)

        subset = df.dropna(subset=[self.setter["selector"]])
        
        if reorder:
            order_constributions = self.reorder_talks(subset)

        else:
            order_constributions = []
            for index, contri in subset.iterrows():                
                order_constributions.append(contri)

        # start text
        text = self.section_separation(section)

        text += "\\input{toc}\n\\newpage\n" if self.toc else ""
        for index, contri in enumerate(order_constributions):
            text += self.contribution(contri, section, index + 1)
        
        with open(output, 'w') as outfile:
            outfile.write(text)
        
        if toc:
            self.poster_TOC('toc.tex')

        
        return text

    def poster_TOC(self, output):
        text = """
\\begin{longtable}{c m{10.5cm} c}
    \\color{hitsblue}{ \\textbf{Stand} } & 
    \\centering \\color{hitsblue}{\\textbf{Poster}} & 
    \\color{hitsblue}{\\textbf{Page}} \\\\[10pt]
    \\endhead
"""
        for poster in self.order_posters:
            stand, title, author = self.order_posters[poster]
            text += f"    {stand} & ``{title}'', {author} &" + \
                    f" \\pageref{{{poster}}} \\\\[17 pt]\n"
        text += "\\end{longtable}\n\\newpage\n"

        with open(output, 'w') as outfile:
            outfile.write(text)
        
        return text
