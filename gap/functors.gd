# GAP Implementation
# This file was generated from 
# $Id: functors.gd,v 1.5 2012/02/27 12:26:34 sunnyquiver Exp $
DeclareAttribute( "DualOfModule", IsPathAlgebraModule );
DeclareOperation( "TransposeOfModule", [IsPathAlgebraModule ] );
DeclareOperation( "DTr", [IsPathAlgebraModule ] );
DeclareOperation( "TrD", [IsPathAlgebraModule ] );
DeclareSynonym("DualOfTranspose", DTr);
DeclareSynonym("TransposeOfDual", TrD);
