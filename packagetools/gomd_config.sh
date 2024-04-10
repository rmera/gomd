#This has been tested in Bash, though I'm rather confident that it will go on zsh as well.


#Not strictly needed, but nice to have
export PATH=$GOMDROOT:$PATH
export PATH=$GOMDROOD/plot:$PATH


#xdrlib
export LD_LIBRARY_PATH=$GOMDROOT/xdrlib/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$GOMDROOT/xdrlib/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$GOMDROOT/xdrlib/include:$CPLUS_INCLUDE_PATH
