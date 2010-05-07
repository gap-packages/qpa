###############################################################################
# HOPF Project Makefile
# DESCRIPTION: This makefile bootstraps the building process of the Hopf project
#
# Copyright, 1998 Virginia Polytechnic Institute and State University.
# Copyright, 1998 Virginia Tech Hopf Project. All rights reserved.
#
# This file may be distributed in accordance with the stipulations existing
# in the LICENSE file accompanying this software.
#
# $Id: Makefile.boot,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
###############################################################################
NOWEBSRCS=hopf.nw legal.nw makefile.nw

Makefile: makefile.nw legal.nw
	-mv -f Makefile Makefile.old
	notangle -t8 -R"[[Makefile]]" $(NOWEBSRCS) > Makefile.tmp
	unexpand Makefile.tmp > Makefile
	rm -f Makefile.tmp
