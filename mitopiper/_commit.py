# The commit this copy of mitopiper was made from. `git archive` (the mitopiper-latest.tar.gz release)
# fills these in (export-subst in .gitattributes); in a git checkout they stay as they are and
# mitopiper.build_info() asks git instead.
COMMIT = "$Format:%h$"
DATE = "$Format:%cs$"
