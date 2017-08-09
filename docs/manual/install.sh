#!/bin/bash

#root directory of project and config file for doxygen
parentname="viennats-dev"
configfile="ViennaTS"

#filename for html shortcut to Documentation outside of html folder
docmainfile="ViennaTS_Doc.html"


#ViennaTS directory is automatically found by searching current working directory
viennaTSdir=$(pwd | sed 's/\(.*viennats-dev\).*/\1/' | sed -e 's/[\/&]/\\&/g')

#list of all commands and arguments that need to be changed for different computers
commands[0]="USE_MDFILE_AS_MAINPAGE"
commands[1]="INPUT"
commands[2]="EXCLUDE"
commands[3]="PROJECT_LOGO"

argument[0]="$viennaTSdir\/README.md"
argument[1]="$viennaTSdir" 
argument[2]="$viennaTSdir\/build $viennaTSdir\/tools $viennaTSdir\/examples"
argument[3]="$viennaTSdir\/docs\/manual\/logo_px55.png"


for i in `seq 0 $(( ${#commands[*]} - 1 ))`; do
	command='sed -i "s/\(.*'${commands[$i]}' *=\).*/\1 '${argument[$i]}'/" '$configfile
	eval $command
done

#run doxygen to make the page & make shortcut outside
doxygen ViennaTS

if [ -e $docmainfile ]
then
	echo "$docmainfile already exists: Not creating."
else
echo '<!DOCTYPE HTML>
<html lang="en-US">
    <head>
        <meta charset="UTF-8">
        <meta http-equiv="refresh" content="1; url=html/index.html">
        <script type="text/javascript">
            window.location.href = "html/index.html"
        </script>
        <title>Page Redirection</title>
    </head>
    <body>
        If you are not redirected automatically, follow this <a href="html/index.html">link to example</a>.
    </body>
</html>' > $docmainfile
fi

