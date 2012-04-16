# GAP Implementation
# This file was generated from 
# $Id: functors.gd,v 1.6 2012/04/16 09:20:23 sunnyquiver Exp $
DeclareAttribute( "DualOfModule", IsPathAlgebraMatModule );
DeclareOperation( "TransposeOfModule", [IsPathAlgebraMatModule ] );
DeclareOperation( "DTr", [IsPathAlgebraMatModule ] );
DeclareOperation( "TrD", [IsPathAlgebraMatModule ] );
DeclareSynonym("DualOfTranspose", DTr);
DeclareSynonym("TransposeOfDual", TrD);
