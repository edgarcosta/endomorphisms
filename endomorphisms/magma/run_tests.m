// usage: magma target:=SUBSTRING run_tests.m
if assigned filename then
  tests := [filename];
else
  tests := Split(Pipe("ls Tests", ""), "\n");
end if;
AttachSpec("spec");
failed := [];
if not assigned target then
  target := "";
end if;

for filename in tests do
  if target in filename then
    fullPath := "tests/" cat filename;
    timestamp := Time();
    try
      printf "%o: ", filename;
      assert eval (Read(fullPath) cat  "return true;");
      printf "Success! %o s\n", Time(timestamp);
    catch e
      Append(~failed, filename);
      printf "Fail! %o s\n %o\n", e, Time(timestamp);;
    end try;
  end if;
end for;
if #failed gt 0 then
  print "Tests failed:";
  for f in failed do
    print f;
  end for;
end if;
if assigned exitsignal then
  exit #failed;
end if;
