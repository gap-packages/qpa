SetPackageInfo( rec(
PackageName := "QPA",
Subtitle := "Quivers and Path Algebras",
Version := "1.35",
Date := "04/01/2024", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

ArchiveURL := Concatenation( "https://folk.ntnu.no/oyvinso/QPA/qpa-",~.Version),

ArchiveFormats := ".tar.gz",

Persons := [
  rec( 
    LastName      := "Green",
    FirstNames    := "Edward",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "green@math.vt.edu",
    WWWHome       := "http://www.math.vt.edu/people/green",
    PostalAddress := Concatenation( [
		       "Department of Mathematics\n",
		       "Virginia Polytechnic Institute and State  University\n",
		       "Blacksburg, Virginia\n",
                       "U.S.A." ] ),
    Place         := "Blacksburg",
    Institution   := "Virginia Polytechnic Institute and State  University"
           ),
  rec( 
    LastName      := "Solberg",
    FirstNames    := "Oeyvind",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "oyvind.solberg@ntnu.no",
    WWWHome       := "https://folk.ntnu.no/oyvinso/",
    PostalAddress := Concatenation( [
		       "Department of Mathematical Sciences\n",
		       "NTNU\n",
		       "N-7491 Trondheim\n",
                       "Norway" ] ),
    Place         := "Trondheim",
    Institution   := "Norwegian University of Science and Technology"
  )              
],

Status := "deposited",

##  You must provide the next two entries if and only if the status is 
##  "accepted" because is was successfully refereed:
# format: 'name (place)'
# CommunicatedBy := "Mike Atkinson (St. Andrews)",
#CommunicatedBy := "",
# format: mm/yyyy
# AcceptDate := "08/1999",
#AcceptDate := "",

README_URL := "https://folk.ntnu.no/oyvinso/QPA/README",
PackageInfoURL := "https://folk.ntnu.no/oyvinso/QPA/PackageInfo.g",

SourceRepository := rec( 
  Type := "git", 
  URL := "https://github.com/gap-packages/qpa"
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),

AbstractHTML := "The <span class=\"pkgname\">QPA</span> package provides data structures \
                   and algorithms for doing computations with finite dimensional quotients \
                   of path algebras, and finitely generated modules over such algebras. The \
                   current version of the QPA package has data structures for quivers, \
                   quotients of path algebras, and modules, homomorphisms and complexes of \
                   modules over quotients of path algebras.",
                   
PackageWWWHome := "https://folk.ntnu.no/oyvinso/QPA/",
               
PackageDoc := rec(
  BookName  := "QPA",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Quivers and Path Algebras",
  Autoload  := true
),

Dependencies := rec(
  GAP := ">=4.5",
  NeededOtherPackages := [["GBNP", ">=0.9.5"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
                      
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := ["quiver","path algebra"],
));
