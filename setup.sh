#!/bin/sh
#
# Run this script from within the root directory of the git clone of your package

set -e

# TODO: error early on if there already is a gh-pages dir
# TODO: error early on if there already is a gh-pages branch
# TODO: make it configurable how/which gap is used
# TODO: make it configurable which remote is used (default: origin)
# TODO: add /gh-pages/ to .gitignore if it is not already in there

# TODO: set UseWorktree based on git version
UseWorktree=No

if [[ ${UseWorktree} = Yes ]]
then
   # Add a new remote pointing to the GitHubPagesForGAP repository
   git remote add -f gh-gap https://github.com/gap-system/GitHubPagesForGAP || :

   # Create a fresh gh-pages branch from the new remote
   git branch gh-pages gh-gap/gh-pages --no-track

   # Create a new worktree and change into it
   git worktree add gh-pages gh-pages
    cd gh-pages
else
   # Create a fresh clone of your repository, and change into it
   url=$(git remote get-url origin)
   git clone ${url} gh-pages
   cd gh-pages

   # Add a new remote pointing to the GitHubPagesForGAP repository
   git remote add gh-gap https://github.com/gap-system/GitHubPagesForGAP
   git fetch gh-gap

   # Create a fresh gh-pages branch from the new remote
   git checkout -b gh-pages gh-gap/gh-pages --no-track
fi

cp -f ../PackageInfo.g ../README* .

git rm -rf doc
mkdir -p doc/
cp -f ../doc/*.{css,html,js,txt} doc/ || :

if [[ -d ../htm ]]
then
  cp -r ../htm .
fi

gap update.g

git add .
# git commit -m "Setup gh-pages based on GitHubPagesForGAP"
# git push --set-upstream origin gh-pages
