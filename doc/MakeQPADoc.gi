# This is to produce the three formats for the QPA-documentation,
# HTML, TXT and PDF formats. The GAPDoc-package is needed to use
# this. 
path := ".";
main := "qpadocumentation.xml";
files:= [];
bookname := "QPA";
MakeGAPDocDoc(path,main,files,bookname);
