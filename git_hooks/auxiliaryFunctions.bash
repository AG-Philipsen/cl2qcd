# This code has been taken from BaHaMAS and adapted
# by Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
#
# Copyright (c) 2017 Alessandro Sciarra
#
# This file is part of CL2QCD.
#
# CL2QCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CL2QCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------#

function errecho() {
    local indentation="   "
    if [ $# -eq 1 ]; then
        printf "$indentation$1\e[0m" 1>&2
    elif [ $# -eq 2 ]; then
        printf "\e[38;5;$2m$indentation$1\e[0m" 1>&2
    elif [ $# -eq 3 ]; then
        printf "\e[$2;38;5;$3m$indentation$1\e[0m" 1>&2
    fi
}

function PrintHookFailure() {
    errecho '\n'
    errecho "HOOK FAILURE ($(basename $0)):" 1 9
    errecho "$@\n" 9
    errecho '\n'
}

function AbortCommit() {
    PrintHookFailure "$1"; shift
    [ $# -gt 0 ] && "$@"
    exit $CL2QCD_failureExitCode
}

#------------------------------------#
# commit-msg hook specific functions #
#------------------------------------#

function GiveAdviceToResumeCommit() {
    errecho 'To resume editing your commit message, run the command:\n\n' 202
    errecho '   git commit -e -F '"$commitMessageFile\n\n" 11 #Variable commitMessageFile from invoking script
}

function IsCommitMessageEmpty() {
    [ -s "$1" ] && return 1 || return 0
}

function RemoveTrailingSpacesAtBeginOfFirstThreeLines() {
    sed -i '1,3{s/^[[:blank:]]*//}' "$1"
}

function RemoveTrailingSpacesAtEndOfEachLine() {
    sed -i 's/[[:blank:]]*$//g' "$1"
}

function AddEndOfLineAtEndOfFileIfMissing() {
    sed -i '$a\' "$1"
}

function CapitalizeFirstLetterFirstLine() {
    sed -i '1s/^\(.\)/\U\1/' "$1"
}

function RemovePointAtTheEndFirstLine() {
    sed -i '1s/[.!?]\+$//g' "$1"
}

function IsFirstLineNotStartingWithLetter() { #Assume no trailing spaces, since removed
    [ $(cat "$1" | head -1 | grep -c '^[[:alpha:]]') -gt 0 ] && return 1 || return 0
}

function IsFirstLineTooShort() {
    [ $(cat "$1" | head -1 | grep -c '^.\{7\}') -gt 0 ] && return 1 || return 0
}

function IsFirstLineTooLong() {
    [ $(cat "$1" | head -1 | grep -c '^..\{60\}') -gt 0 ]  && return 0 || return 1
}

function IsSecondLineNotEmpty() {
    [ $(wc -l < "$1") -lt 2 ] && return 1 #Needed otherwise head and tail below match first line
    [ $(cat "$1" | head -2 | tail -1 | grep -c '^[[:blank:]]*$') -gt 0 ]  && return 1 || return 0
}

function IsAnyOfTheLinesAfterTheSecondTooLong() {
    [ $(cat "$1" | tail -n +2 | grep -c '^..\{72\}') -gt 0 ] && return 0 || return 1
}

#------------------------------------#
# pre-commit hook specific functions #
#------------------------------------#

function GiveAdviceAboutUserNameAndEmail() {
    errecho 'Use the commands\n' 202
    errecho '   git config --global user.name "Your Name"\n' 11
    errecho '   git config --global user.email "you@yourdomain.com"\n' 11
    errecho 'to introduce yourself to Git before committing.\n\n' 202
    errecho 'Omit the "--global" option to set your infortmation only in the local repository.\n\n' 208
}

function GiveAdviceAboutUserNameFormat() {
    errecho 'Please, configure your user.name using the command\n' 202
    errecho '   git config --global user.name "Your Name"\n' 11
    errecho 'where your name has to be formed by two word starting with\n' 202
    errecho 'capital letter and separated by one space.\n\n' 202
    errecho 'Omit the "--global" option to set your infortmation only in the local repository.\n\n' 208
}

function GiveAdviceAboutUserEmailFormat() {
    errecho 'Please, configure your user.email using the command\n' 202
    errecho '   git config --global user.email "you@yourdomain.com"\n' 11
    errecho 'where your email has to be in a valid format as shown here above.\n\n' 202
    errecho 'Omit the "--global" option to set your infortmation only in the local repository.\n\n' 208
}

function GiveAdviceAboutNonASCIICharacters() {
    errecho 'This can cause problems if you want to work with people on other platforms.\n' 202
    errecho 'To be portable it is advisable to rename the file.\n' 202
    errecho 'If you know what you are doing you can disable this check using:\n' 202
    errecho '   git config hooks.allownonascii true\n\n' 11
}

function GiveAdviceAboutBranch() {
    errecho "If you are sure about what you are doing, commit with \e[38;5;11m--no-verify\e[38;5;202m to bypass this check.\n\n" 202
}

function GiveAdviceAboutWhitespaceError() {
    errecho 'Use the command\n' 202
    errecho "   git diff-index --check --cached $againstSHAToCompareWidth\n" 11
    errecho 'to have a look to the whitespace violation on staged files.\n\n' 202
}

function GiveAdviceAboutAstyle() {
    errecho 'Please install it before continuing\n' 202
    errecho '  http://astyle.sourceforge.net/\n\n' 11
}

function PrintReportOnFilesWithStyleErrors() {
    errecho "=================================================================================================\n" 14
    errecho "   Here a list of the affected files:\n" 202
    for file in "$@"; do
        errecho "     - ${file}\n" 208
    done
    errecho "\n"
    errecho " Please fix before committing. Don't forget to run git add before trying to commit again.\n" 202
    errecho " If the whole file is to be committed, this should work (run from the top-level directory):\n" 202
    errecho "\n"
    errecho "    astyle ${astyleParameters} $file; git add $file; git commit" 11
    errecho "\n"
    errecho "=================================================================================================\n" 14
}
