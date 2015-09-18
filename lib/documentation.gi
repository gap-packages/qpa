InstallGlobalFunction( MakeQPADocumentation,
function()
  local dir, path, main, files, bookname;
  dir := DirectoriesPackageLibrary( "QPA", "doc" );
  path := Filename( dir[1], "" );
  main := "qpadocumentation.xml";
  files:= [];
  bookname := "QPA";
  MakeGAPDocDoc( path, main, files, bookname, "MathJax" );
  CopyHTMLStyleFiles( path );
end );
