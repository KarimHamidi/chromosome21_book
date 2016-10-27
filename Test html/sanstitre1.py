# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 14:03:32 2016

@author: karim
"""

import os 
os.chdir('C:\\Users\\karim\\FirstStepProject\\.spyproject\\Test html')

h = open('test.html','w')
message = """<html>
    <head>
   <title>My first HTML document</title>
    </head>
    <body>
   <p>Hello <strong>world!</strong><p>
   <br>
   <p><embed src="Seqmin80.txt" width="2000" height="2000"><p>
    </body>
</html>"""
h.write(message)
h.close()
import pdfkit
pdfkit.from_url('file:///C:/Users/karim/FirstStepProject/.spyproject/Test%20html/test.html', 'micro.pdf')
import pdfkit

config = pdfkit.configuration(wkhtmltopdf=path_wkthmltopdf)
pdfkit.from_url("http://google.com", "out.pdf", configuration=config)
import nbconvert
import sys
from bs4 import BeautifulSoup
