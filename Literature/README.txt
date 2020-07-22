###############################################################################
READ THE FOLLOWING DOCUMENT CAREFULLY BEFORE MAKING CHANGES TO THIS FOLDER!!!
###############################################################################

# ---------------------------
# Header
# ---------------------------
  
Peridynamic literature README

Contact:   Martin Rädel,  martin.raedel@dlr.de
           DLR Composite Structures and Adaptive Systems
           Mail:  martin.raedel@dlr.de
           Phone: +49 (0) 531 295-2048

                                __/|__
                               /_/_/_/  
           www.dlr.de/fa/en      |/ DLR

# ---------------------------
# Concept
# ---------------------------

- The .bib-file shall contain any document related to peridynamics
- It can be loaded to any document
- If you need other references to non-peridynamic literature, please create an individual .bib-file in your project
- Do not add any non-peridynamic related publications in this .bib-file

# ---------------------------
# The .bib-file
# ---------------------------

- Several projects use the file "PeridynamicLiterature.bib" from this folder
- Do's:
    - In case you add a pdf to one of the Literature subfolders, please add the proper BibTeX-entry in the .bib-file right away
    - Usually literature sources allow some kind of BibTeX-export:
        - Use it if available
        - Adapt the formatting to the existing entries
        - Add the document:
            - In the proper BibTeX-subcategory
            - Sorted alphabetically according to BibTeX-Key
       - Check if the:
            - Document information is as complete as it gets
            - BibTeX-key for the new document is not used already
    - Test ALL projects for proper compilation after you made changes to the document
- Dont's:
    - Do not mess around with this file unless you know what you are doing!
    - NEVER EVER EVER change existing BibTeX-keys
- Conventions:
    - BibTeX-Key:
        1. Last name of the first author
        2. Initials of first author first name
        3. Year of publication
        4. Optional: Counter a,b,c... in case of multiple identical combinations of 1.-3.
    - If known, type author first names in full length in Author entry
    - Do not save the abstract
    - Add a comma after each entry, even the last one, for copy-paste-reasons
    - Convert any umlaut or other encoding-dependent symbols to the corresponding base LaTeX command (ä -> \"a, ö -> \"o, ü -> \"u, ß -> \ss{})
