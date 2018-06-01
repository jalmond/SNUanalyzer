export SNUVERSION=v9-4-1
### If there is a small bug/new code then new subtag is made
export tag_numerator='.1'
if [[ $1 == 'branch' ]];
    then
    export SNUTAG=
else
    export SNUTAG=v9-4-1.1
fi

