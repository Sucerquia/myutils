import pandas as pd
from myutils.miscellaneous import output_terminal


def add_day(day):
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

def session_title(event):
    """
    Add the title of the session including the chair before showing the events
    of the session.

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

def add_item(event, timespace=1, titlespace=10.5):
    time = event["startTime"][17:22]
    title = event["title"].replace("&#8211;", "--").replace("\\/", "/")
    text = f"""
\\begin{{tabular}}{{p{{{timespace}cm}}p{{{titlespace}cm}}}}
    {time} & {title} \\\\
\\end{{tabular}}
"""
    return text

def add_program():
    """
    Copy the schedule from the webpage and create the program.

    Return
    ======
    (str) Program in LaTeX format.
    """
    all_events = output_terminal('curl https://simplaix-workshop2025.h-its.org ' +\
                                '| grep \'"items"\' | sed "s/<\/script>//g" ',
                                print_output=False)
    all_events = eval(all_events)

    events_per_day = {"Wednesday, 7 May 2025": [],
                    "Thursday, 8 May 2025": [],
                    "Friday, 9 May 2025": []}

    for event in all_events["items"]:
        if "Wed, 07 May 2025" in event["startTime"]:
            events_per_day["Wednesday, 7 May 2025"].append(event)
        elif "Thu, 08 May 2025" in event["startTime"]:
            events_per_day["Thursday, 8 May 2025"].append(event)
        elif "Fri, 09 May 2025" in event["startTime"]:
            events_per_day["Friday, 9 May 2025"].append(event)
    text=""

    for day in events_per_day:
        text += add_day(day)
        session = False
        for event in events_per_day[day]:
            if "Session" ==  event["title"][:7]:
                if "&#8211; Chair:" in  event["title"]:
                    text += session_title(event)
                    session=True
                else:
                    title_and_kind = event["title"].split("\\/")[1]
                    kind = title_and_kind.split(":")[0]
                    author = title_and_kind.replace(kind + ": ", "")
                    if event['excerpt'] != '':
                        title = event['excerpt'].replace("&#8220;",
                                                        "``").replace("&#8221;",
                                                                    "''")
                    event["title"] = f" ({kind}) {author} -- {title}"
                    text += add_item(event)
                    continue
            else:
                if session:
                    text += "\n\\vspace{0.25 cm}\n"
                    session = False
                text += "\n\\vspace{0.25 cm}\n"
                text += add_item(event)
        text += "\n}\n"
    
    return text

def section_separation(title):
    """
    Adds a new page separating the sections.

    Parameters
    ==========
    title: str
        title of the section.
    
    Return
    ======
    (str) A new page separating sections with a big title in an isolated page.
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

def title_contri(title):
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
    \\textbf{{\\Large \\color{{hitsblue}}{{ {costume_title(title)} }} }}
\\end{{center}}

"""
    return text

def author_line(authors_info, affiliation):
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

def costume_title(text):
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
    words that should not be transformed such as acronyms. Access and modify
    according to your preferences. 
    """
    exceptions = {"a", "an", "and", "as", "at", "but", "by", "for", "in",
                  "nor", "of", "on", "or", "so", "the", "to", "up", "yet",
                  "with", "through", "into"}
    acronyms = {"ART-SM:", "ML-based", "MACE", "ML/MM", "QM/MM", "OF-DFT",
                "LHCII", "AI"}
    words = text.split()
    titled = [words[0] if words[0] in acronyms else
              words[0].capitalize()]  # Always capitalize the first word
    for word in words[1:]:
        word = word if ((word in exceptions) or (word in acronyms)) \
                        else word.capitalize()
        # capitalize also after -
        if '-' in word and word not in acronyms:
            start = [word[0]]
            word = ''.join(start + [word[i].upper() if word[i - 1] == '-'
                                    else word[i] for i in range(1, len(word))])
        titled.append(word)
    return ' '.join(titled)


def contribution(row, kind, stand=None):
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
    if kind == 'talk':
        text = title_contri(costume_title(row["Oral title"]))
    elif kind == 'poster':
        text = title_contri(costume_title(row["Poster title"]))
    else:
        raise ValueError("Wrong type of contribution")

    name = row["firstName"].title()
    lastname = row["lastName"].title()
    email = row["email"]
    complete_name = name + " " + lastname + " (" + email + ')'
    if row['title'] == 'Dr.' or row['title'] == 'Dr':
        complete_name = 'Dr. ' + complete_name
    affiliation = row["institute"] + ', ' + row["city"] + ', ' + row['country']
    text += author_line(complete_name, affiliation)
            

    if kind == 'talk':
        text += row["oral abstract"] + '\n\n'
    if kind == 'poster':
        text += row["Poster abstract"] + '\n\n'
        text += '\\begin{flushright}\n' + \
                f'    Stand \\Ntag {stand}\n' + \
                '\\end{flushright}\n'

    text += '\\newpage\n\n'
    return text


# add2executable
def abstract_book(excel_file, output):
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

    text = """
\\documentclass{article}
\\usepackage{xcolor}
\\usepackage{setspace}
\\usepackage{graphicx} % Required for inserting images

\\definecolor{hitsblue}{HTML}{1f4a87}

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

"""
    text += section_separation('Program')
    text += add_program()

    df = pd.read_excel(excel_file)
    # Display the first few rows
    posters = df.dropna(subset=['Poster'])
    talks = df.dropna(subset=['Oral contribution'])

    text += section_separation('Talks')
    for index, contri in talks.iterrows():
        text += contribution(contri, 'talk')
    
    text += section_separation('Posters')
    
    for index, contri in posters.iterrows():
        text += contribution(contri, 'poster', stand=index + 1)
    
    text += '\n\\end{document}\n'
    
    with open(output, 'w') as outfile:
        outfile.write(text)
