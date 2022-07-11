LoadPackage("qpa");
dir := DirectoriesPackageLibrary("qpa", "tst");

# compareFunction: performing "uptowhitespace", and ignoring the incorrect
# pluralisation "1 generators".
# See https://github.com/gap-packages/qpa/issues/71
func := function(a, b)
  # Remove newlines
  a := ReplacedString(ShallowCopy(a), "\\\n", "");
  b := ReplacedString(ShallowCopy(b), "\\\n", "");
  # Remove whitespace
  RemoveCharacters(a, " \n\t\r");
  RemoveCharacters(b, " \n\t\r");
  # Remove (probably) incorrect pluralisation
  a := ReplacedString(a, "1generators", "1generator");
  b := ReplacedString(b, "1generators", "1generator");
  return a = b;
end;

TestDirectory(dir, rec(exitGAP := true,
                       testOptions:=rec(compareFunction:=func)));

FORCE_QUIT_GAP(1);
