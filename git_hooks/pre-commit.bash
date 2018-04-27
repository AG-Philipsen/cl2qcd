#!/bin/bash
#
# Check that the code follows a consistent code style and that
# some given rules are respected by what is about to be committed.
#
# LIST OF CHECKS:
#  1) Check if user name and email are set and reasonable
#  2) Check branch on which the commit is being done
#  3) Prevent non-ASCII characters in filenames
#  4) Make whitespace git check
#  5) Go through fully staged files and
#      - remove trailing spaces at end of lines
#      - add end of line at the end of file
#      - remove empty lines at the end of the file
#  6) Go through staged files and
#      - check copyright statement
#      - check license notice
#      - check code style using clang-format
#
# Called by "git commit" with no arguments.  The hook should
# exit with non-zero status after issuing an appropriate message if
# it wants to stop the commit.
#
# The code style check  has been "stolen" from PyGit and slightly
# modified by Matthias Bach <bach@compeng.uni-frankfurt.de>
#
# The rest of the hook has been taken from BaHaMAS and adapted
# by Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
#
# Copy this script to .git/hooks to enforce proper coding style on
# committed files (or use automatic hooks setup)
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2017,2018 Alessandro Sciarra
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

readonly CL2QCD_hooksFolder="$(dirname "$(readlink -f "$BASH_SOURCE")")"
readonly CL2QCD_failureExitCode=1
source ${CL2QCD_hooksFolder}/auxiliaryFunctions.bash || exit $CL2QCD_failureExitCode

#-----------------------------------------------------------------------------
#Check for committer identity as such and check if it exists in history
if [ "$GIT_AUTHOR_NAME" != '' ]; then
    readonly userName="$GIT_AUTHOR_NAME"
else
    readonly userName="$(git config --get user.name)"
fi
if [ "$GIT_AUTHOR_EMAIL" != '' ]; then
    readonly userEmail="$GIT_AUTHOR_EMAIL"
else
    readonly userEmail="$(git config --get user.email)"
fi
if [ "$GIT_COMMITTER_NAME" != '' ]; then
    readonly committerName="$GIT_COMMITTER_NAME"
else
    readonly committerName="$(git config --get user.name)"
fi
if [ "$GIT_COMMITTER_EMAIL" != '' ]; then
    readonly committerEmail="$GIT_COMMITTER_EMAIL"
else
    readonly committerEmail="$(git config --get user.email)"
fi

if [ "$userName" = '' ] || [ "$userEmail" = '' ]; then
    AbortCommit "User information not configured!" GiveAdviceAboutUserNameAndEmail
fi
if [[ ! $userName =~ ^[A-Z][a-z]+(\ [A-Z][a-z]+)+$ ]]; then
    AbortCommit "User name not allowed." GiveAdviceAboutUserNameFormat
fi
if [[ ! $userEmail =~ ^[^@]+@[^@]+$ ]]; then
    AbortCommit "User email not allowed." GiveAdviceAboutUserEmailFormat
fi
if [[ ! $committerName =~ ^[A-Z][a-z]+(\ [A-Z][a-z]+)+$ ]]; then
    AbortCommit "Committer name not allowed." GiveAdviceAboutCommitterNameFormat
fi
if [[ ! $committerEmail =~ ^[^@]+@[^@]+$ ]]; then
    AbortCommit "Committer email not allowed." GiveAdviceAboutCommitterEmailFormat
fi

readarray -t existingUserOrCommitterNamesAndEmails <<< "$(cat <(git log --all --format='%an %ae') <(git log --all --format='%cn %ce') | sort -u)"
readonly existingUserOrCommitterNamesAndEmails
newUser=0
for userOrCommitterNameAndEmail in "${existingUserOrCommitterNamesAndEmails[@]}"; do
    if [ "$userOrCommitterNameAndEmail" = "$userName $userEmail" ] || [ "$userOrCommitterNameAndEmail" = "$committerName $committerEmail"]; then
        newUser=1
        break
    fi
done
if [ $newUser -eq 0 ]; then
    errecho "\n"
    errecho '\e[1mWARNING:\e[21m A new author and/or committer name(s) and/or email(s) has(have) been found.\n' 202
    errecho "              Author:   $userName  <$userEmail>\n" 11
    errecho "           Committer:   $committerName  <$committerEmail>\n\n" 11
    AskUser "Would you like to continue the commit creating, then, a new author ot committer?" 202
    if UserSaidNo; then
        AbortCommit "Author and/or committer name(s) and/or email(s) to be checked!"
    fi
fi

#-----------------------------------------------------------------------------
#Check branch: direct commit on 'develop' should not be done
readonly actualBranch="$(git rev-parse --abbrev-ref HEAD)"
if [ "$actualBranch" = 'develop' ]; then
    AbortCommit "Direct commits on \"$actualBranch\" branch should be avoided!" GiveAdviceAboutBranch
fi

#-----------------------------------------------------------------------------
#Check added filenames (modified from pre-commit.sample)
readonly headSHA="$(git rev-parse --verify HEAD 2>/dev/null)"
if [ "$headSHA" != '' ]; then
    readonly againstSHAToCompareWidth=HEAD
else
    # Initial commit: diff against an empty tree object
    readonly againstSHAToCompareWidth=$(git hash-object -t tree /dev/null) # it gives 4b825dc642cb6eb9a060e54bf8d69288fbee4904
fi

# If you want to allow non-ASCII filenames set this variable to true.
readonly allownonascii=$(git config --bool hooks.allownonascii)

# Cross platform projects tend to avoid non-ASCII filenames; prevent
# them from being added to the repository. We exploit the fact that the
# printable range starts at the space character and ends with tilde.
#
# NOTE: the use of brackets around a tr range is ok here, (it's
#       even required, for portability to Solaris 10's /usr/bin/tr), since
#       the square bracket bytes happen to fall in the designated range.
if [ "$allownonascii" != 'true' ] && [ $(git diff --cached --name-only --diff-filter=A -z $againstSHAToCompareWidth | LC_ALL=C tr -d '[ -~]\0' | wc -c) -ne 0 ]; then
    AbortCommit 'Attempt to add a non-ASCII file name.' GiveAdviceAboutNonASCIICharacters
fi

#-----------------------------------------------------------------------------
#Do not allow spaces in filenames
if [ $(git diff --cached --name-only --diff-filter=A -z $againstSHAToCompareWidth | grep -c '[ ]') -ne 0 ]; then
    AbortCommit "Spaces are not allowed in filenames!"
fi

#-----------------------------------------------------------------------------
# Get list of staged files which have to be checked with respect to CODE STYLE
readonly listOfFilesToBeChecked=( $(git diff-index --cached --name-only HEAD --diff-filter=ACMR | grep -e "\.c$" -e "\.cpp$" -e "\.h$" -e "\.hpp$" -e "\.cl$") )
if [ ${#listOfFilesToBeChecked[@]} -ne 0 ]; then
    # Check for existence of clang-format, and error out if not present
    if ! builtin type -P clang-format >>/dev/null; then
        AbortCommit "The program \"clang-format\" was not found!" GiveAdviceAboutClangFormat
    else
        errecho "\n"
        errecho "Checking style of code... " 39
        readonly clangFormatParameters="-style=file"
        filesWithCodeStyleErrors=()
        for file in "${listOfFilesToBeChecked[@]}"; do
            # The file in the staging area could differ from the ${file} in case
            # it was only partially staged. It would be possible, in principle, to
            # obtain the file from the index using
            #     indexFile=$(git checkout-index --temp ${file} | cut -f 1)
            #
            # Actually, if we use a copy of ${file} with a different name
            # in order to check if the code style rules are respected,
            # clang-format could reorder differently the #include, since
            # it might not identify that the ${file} is a .cpp associated
            # with a .hpp (in this case the .hpp should be included as first).
            # Therefore, we always check the ${file} as if it was fully staged
            # and if there are style violation in the not staged modifications
            # we do not care and we abort the commit any way (just give some
            # advice to the user).
            newfile=$(mktemp /tmp/$(basename ${file}).XXXXXX) || exit ${CL2QCD_failureExitCode}
            clang-format ${clangFormatParameters} ${file} > ${newfile} 2>> /dev/null
            if ! cmp -s "${file}" "${newfile}"; then # <- if they differ
                filesWithCodeStyleErrors+=( "${file}" )
            fi
            rm "${newfile}"
        done
        if [ ${#filesWithCodeStyleErrors[@]} -ne 0 ]; then
            echo; AbortCommit "Code style error found!" PrintReportOnFilesWithStyleErrors "${filesWithCodeStyleErrors[@]}"
        else
            errecho " done!\n" 39
        fi
    fi
fi

#-----------------------------------------------------------------------------
# Work on tab/spaces in files (only those fully stages, since after
# modification we have to add them and if done on partially staged
# they would be then fully staged against user willing!
# (adapted from https://gist.github.com/larsxschneider/3957621)
# NOTE: Assume no spaces in filenames (checked above)
readonly stagedFiles=( $(git diff-index --name-only --cached --diff-filter=AM $againstSHAToCompareWidth | sort | uniq) )
readonly partiallyStagedFiles=( $(git status --porcelain --untracked-files=no | # Find all staged files
                                         egrep -i '^(A|M)M '                  | # Filter only partially staged files
                                         sed -e 's/^[AM]M[[:space:]]*//'      | # Remove leading git info
                                         sort | uniq) )                         # Remove duplicates

# Merge staged files and partially staged files
readonly stagedAndPartiallyStagedFiles=( "${stagedFiles[@]}" "${partiallyStagedFiles[@]}" )

# Remove all files that are staged *AND* partially staged -> we get only the fully staged files
readonly fullyStagedFiles=( $(tr ' ' '\n' <<< "${stagedAndPartiallyStagedFiles[@]}" | sort | uniq -u) )

if [ ${#fullyStagedFiles[@]} -ne 0 ]; then
    errecho "\n"
    errecho "Fixing trailing whitespaces and newline at EOF in fully staged files:\n" 39
    for file in "${fullyStagedFiles[@]}"; do
        errecho "   - $file\n" 87

        # Strip trailing whitespace
        sed -i 's/[[:space:]]*$//' "$file"

        # Add newline to the end of the file
        sed -i '$a\' "$file" # 'a\' appends the following text, which is nothing, in this case!
        # The code "$a\" just says "match the last line of the file, and add nothing to it."
        # But, implicitly, sed adds the newline to every line it processes if it is not already there.

        # Remove empty (w/o spaces) lines at the end of the file (http://sed.sourceforge.net/sed1line.txt)
        sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' "$file"
        # Alternative in awk, but it needs a temporary file
        # awk '/^$/{emptyLines=emptyLines"\n"; next} {printf "%s", emptyLines; emptyLines=""; print}'

        # Stage all changes
        git add "$file"
    done
fi

#-----------------------------------------------------------------------------
# Check CL2QCD header in staged files
if [ ${#stagedFiles[@]} -ne 0 ]; then
    readonly CL2QCD_headerFile="${CL2QCD_hooksFolder}/header.txt"
    readonly numberOfExpectedTextLines=$(sed '/^$/d' "${CL2QCD_headerFile}" | wc -l)
    readonly expectedCopyright='Copyright \(c\) ([2][0-9]{3}[,-]?[ ]?)*'"$(date +%Y) ${userName}"
    filesWithWrongOrMissingHeader=()
    filesWithMissingCopyright=()
    errecho "\n"
    errecho "Checking header of staged files..." 39
    for file in "${stagedFiles[@]}"; do
        if [[ ${file} =~ [.](bash|c|C|cl|cmake|cpp|h|hpp|m|nb|py|sh|tex|txt)$ ]]; then
            # Note that here the "| sort | uniq" avoids to count multiple times lines
            # of the CL2QCD_headerFile which could be by accident repeated in third party code
            numberOfMatchingLines=$(grep -o -f "${CL2QCD_headerFile}" "${file}" | sort | uniq | wc -l)
            if [ $numberOfMatchingLines -ne $numberOfExpectedTextLines ]; then
                filesWithWrongOrMissingHeader+=( "${file}" )
            fi
            if [ $(grep -cE "$expectedCopyright" "${file}") -eq 0 ]; then
                filesWithMissingCopyright+=( "${file}" )
            fi
        fi
    done
    if [ ${#filesWithMissingCopyright[@]} -ne 0 ] || [ ${#filesWithWrongOrMissingHeader[@]} -ne 0 ]; then
        errecho "\n"
        if [ ${#filesWithWrongOrMissingHeader[@]} -ne 0 ]; then
            PrintReportOnFilesWithWrongOrMissingHeader "${filesWithWrongOrMissingHeader[@]}"
        fi
        if [ ${#filesWithMissingCopyright[@]} -ne 0 ]; then
            PrintReportOnFilesWithMissingCopyright "${filesWithMissingCopyright[@]}"
        fi
        AskUser "Would you like to continue the commit without fixing the header of the files?"
        if UserSaidNo; then
            AbortCommit "Files with wrong or missing CL2QCD header found!" PrintSuggestionToFixHeader
        fi
    else
        printf " \e[92mdone!\e[0m\n"
    fi
fi

#-----------------------------------------------------------------------------
#If there are still whitespace errors, print the offending file names and fail.
if ! git diff-index --check --cached $againstSHAToCompareWidth >/dev/null 2>&1; then
    AbortCommit "Whitespace errors present in staged files!" GiveAdviceAboutWhitespaceError
fi
