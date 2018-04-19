How to contribute
=================

If you arrived on this page, it is probable that you are going to contribute to CL2QCD and, then, first of all, thanks for your work.
If instead you are here for another purpose, we hope you will find what you search.

CL2QCD grew over the years and it has become quite a large code.
A uniform style and way of behaving by all the authors is crucial to contribute keeping it clean and ordered.
Deciding to contribute to CL2QCD, you agree to abide by our [code of conduct](CODE_OF_CONDUCT.md), which we encourage you to read, as well as to follow the rules described here in the following.


Quick start
-----------

If you are part of the CL2QCD community, just ignore this section.
Still, read carefully the rest of the page and ask other members in case you have any doubt.

To contribute to the project as external contributor, a natural, quite common and suggested way consists of few steps.
* Fork the repository on GitHub (this requires having a [GitHub account](https://github.com/signup/free)).
* Clone it to your local machine and set up the distributed [git hooks](#using-the-distributed-git-hooks).
* Create a branch from where you want to base your work.
  **Avoid** to work directly on the `master` branch!
* Develop the code as you wish regarding the content, but respect the style you find.
* Commit your changes without violating our [commit rules](#commit-rules).
* Be sure to have tested the new code.
  In any case run **all** the existing tests (running `ctest` in your `build` folder) to check that nothing was accidentally broken.
  Test your code on CPUs as well as on GPUs and, if possible, on different devices.
* Push your work to your branch in the remote, forked repository.
* Submit a [pull request](https://help.github.com/articles/about-pull-requests/), leaving a *good* description of your work.
* Wait a feedback by CL2QCD developers.
  Depending on the complexity of your work, this could take more or less time.
  The better you followed our guidelines, the quicker the reviewing process will go.


Uniforming the style
--------------------

### Commit rules

Around in the web, it is plenty of advice about how writing a good commit message.
Reading here and there, we came up with our standard, which we report in the following.

##### Rules that MUST be respected

1. Limit the subject line to **50** characters, but use at least **8**
1. Capitalise the subject line
1. Do not end the subject line with a period
1. Separate subject from body with a blank line
1. Wrap the body at **72** characters

##### Rules that MAY be followed

* Use the imperative mood in the subject line
* Use the body to explain what and why vs. how
* Use [autolinked references](https://help.github.com/articles/autolinked-references-and-urls/) if relevant

We encourage you to read [more about them](https://chris.beams.io/posts/git-commit/).
It should not be so hard to get used to the above rules.
However, to be sure to fulfil them, you should use the [git hooks](#using-the-distributed-git-hooks) distributed in the repository, which can be automatically set up.


### Editing existing files or creating new ones

CL2QCD is distributed under the terms of the GPL.
In general, try to follow the [instructions of the GPL](http://www.gnu.org/licenses/gpl-howto.en.html).
Here in the following, some guidelines for authors are provided in order to reach a coherent style of copyright and licensing notices.

* If you are contributing for the first time, add you name and contact address (email) to the file [AUTHORS](AUTHORS.md).
* If you create a new file, add copyright and license remarks to a comment at the top of the file (after a possible [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) and short introductory comments).
  Templates are provided below.
  Copyright remarks should not only use the **(c)** symbol, but also the term **Copyright**.
* When editing an existing file, ensure there is a copyright remark with your name.
  If there is none, add one directly below any existing copyright remarks.
  If there is already a copyright remark with your name, ensure it contains the current year.
  Otherwise add the current year.
  The years should form a comma-separated list.
  If more than two subsequent years are given, the list can be merged into a range.

The following sample shows how the top level comment of a `C++` source file should look like.

```cpp
/** @file
 * Cool functionality to show some neat things
 *
 * Copyright (c) 2013 Jane Doe <doe@example.com>
 *
 * [LICENSE notice]
 */

```

For a `Python` source file, a shebang is needed and has to be included.

```python
#!/usr/bin/env python
# coding=utf8
#
# Cool functionality to show of some neat things
#
# Copyright (c) 2012, 2013 Jane Doe <doe@example.com>
# Copyright (c) 2013 Max Mustermann <mustermann@example.com>
#
# [LICENSE notice]

```

Please, include an empty line between the comment and the first line of code.
The [`[LICENSE notice]`](git_hooks/header.txt) has to be fully copied and out-commented according to the language of the source file.
Again, to be sure to have done everything according to our standard, we encourage to use our [git hooks](#using-the-distributed-git-hooks), which check about the license notice and the copyright remark automatically.


### Using the distributed git hooks

CL2QCD contains a set of git hooks which are meant to help the developer to respect the common style and to minimise the risk of violating the community rules.
Executing the [`createHooksSymlink.bash`](git_hooks/createHooksSymlink.bash) will set up the available hooks.
Since [`clang-format`](http://releases.llvm.org/) is used to uniform the code style in the remote repository, you need to have it installed.
It is part of the LLVM suite, which can be download as pre-built binaries for your operative system.
Once extracted the `.tar.xz`, make sure that the `bin/clang-format` executable can be automatically found, e.g. setting the `PATH` environment variable accordingly.

In principle, you can use your favourite style developing the code, but you should run `clang-format` on the modified files before committing them.
However, we encourage you to use our style straight away, since most of the editors and IDE offer the possibility to integrate `clang-format`.
Have a look [here](https://clang.llvm.org/docs/ClangFormat.html#vim-integration) for more information, or check out on the web if your editor is not here mentioned.
For example, if you use `Eclipse`, you can take advantage of the [`CppStyle` plug-in](https://github.com/wangzw/CppStyle).

After having run the `createHooksSymlink.bash` script, you will maybe notice that a `_clang-format` symlink will appear at the top-level of the CL2QCD repository.
This is meant to contain the code style options and it has to be placed there in order to be naturally found by `clang-format` (which should be run with the `-style=file` option).

#### `pre-commit`

* It checks that git user name and the git email are set and reasonable.
   * The user name should contain only capitalised words with alphabetic characters, i.e. [[A-Za-z]](https://en.wikipedia.org/wiki/Regular_expression#Character_classes).
   * The email should be a valid one, preferably an institutional one.
* It checks the user and committer identities and ask for confirmation in case they are not present in history.
* It checks the branch name on which the commit is being done.
* It prevents non-ASCII characters in filenames.
* It makes the whitespace git check.
* It goes through fully staged files and
   * it removes trailing spaces at end of lines;
   * it adds an end of line at the end of the file;
   * it removes empty lines at the end of the file.
* It goes through all staged files and
   * it checks the copyright statement;
   * it checks the license notice;
   * it checks the code style using `clang-format` (which must be installed).

#### `commit-msg`

This hooks enforces the [commit rules](#commit-rules), together with some more style check.

* The commit is aborted if
  * the first line of the commit has more than **50** characters;
  * the first line of the commit has less than **8** characters;
  * the first line does not begin with a character;
  * the second line is not empty;
  * any following line after the second has more than 72 characters.
* The commit message is changed according to the following cases,
  * a small letter at the beginning of the first line is capitalised;
  * any character among `.?!` at the end of the first line is repetitively removed (e.g. `commit..!!!` becomes `commit`);
  * trailing spaces at begin and end of the first two lines are removed.

#### Useful information

Most of the editors have the possibility to help you in typing a commit message.

---

**`vim`** does it in a natural way and even with default configuration provides a nice syntax highlighting for commit messages.
Using colours, it notifies you if your commit summary is too long or if you are typing something on the second line that should be empty.
To automatically wrap the commit description at 72 characters, you could add `set textwidth=72` to your `.vimrc` file in your home directory.
To make git use `vim` as editor for all your repositories, use something like `git config --global core.editor vim`.

---

If your favourite editor is **`emacs`**, the way is not so down-hill.
But it is not so tough neither.
From MELPA you can download the [git-commit package](https://melpa.org/#/git-commit), which will provide you with useful functionality.
Then you can add few lines to your `.emacs`, which could look like below.

```lisp
;; Git specific operation to simplify environment when committing
(setq column-number-mode t)              ; show column number in the mode line
(global-git-commit-mode)                 ; activate git-commit mode
(setq git-commit-summary-max-length 49)  ;  - changing style in first line after the 50th char
(setq git-commit-fill-column 71 )        ;  - activate auto-fill-mode on space and return
(aset auto-fill-chars ?. t)              ;    after the 72th char. It can be useful to have
(aset auto-fill-chars ?? t)              ;    line folded also on any of the punctuation like
(aset auto-fill-chars ?! t)              ;    .?! which could be inserted beyond char 72.
```

If your `.emacs` file takes some time to be loaded and you would like to quickly fire up emacs to make a commit, then you could put the code above in some lighter standalone file (e.g. `${HOME}/.emacs_for_git`) and then tell git how to invoke properly the editor (e.g. `git config --global core.editor "emacs -nw -Q -l ${HOME}/.emacs_for_git"`).

---
