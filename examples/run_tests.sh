for file in Test-*; do echo $file; time magma Utils/SetQuitOnErrortrue.m $file Utils/Exit.m> /dev/null|| break; done
