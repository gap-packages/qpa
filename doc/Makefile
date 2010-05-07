# -*- mode: Noweb; noweb-code-mode: makefile-mode -*-%%
###############################################################################
# QPA User Manual Makefile
# DESCRIPTION: Makefile
#
# Copyright, 1998, 1999 Virginia Polytechnic Institute and State University.
# Copyright, 1998, 1999 Virginia Tech Path Algebra Project. All rights reserved.
#
# This file may be distributed in accordance with the stipulations existing
# in the LICENSE file accompanying this software.
#
# $Id: Makefile,v 1.1 2010/05/07 13:16:24 sunnyquiver Exp $
###############################################################################
TEXSRCS = manual.tex

all: clean manual

manual: manual.ps

manual.ps: manual.dvi
	dvips -t letter -o manual.ps manual.dvi

manual.dvi: $(TEXSRCS)
	-tex '\batchmode\input manual' 
#	bibtex manual
#	tex '\batchmode\input manual'
	tex '\scrollmode\input manual'

clean:
	rm -f *.ps *.dvi *.aux *.log *.toc *.bbl *.blg
	rm -f *.lab *.idx *.six
