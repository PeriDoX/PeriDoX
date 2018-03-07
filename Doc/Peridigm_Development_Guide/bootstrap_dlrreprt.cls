%
% This is file `bootstrap_dlrreprt.cls', 
% a bootstrap version of 'dlrreprt.cls'
%
% REVISIONS:    2018-01-31 init     Martin Raedel <martin.raedel@dlr.de>
%                                   Based on dlrreprt-class and styles
%                                   Secondpage information not included
%                                   Can be built and input individually
%               
% Contact       Martin Raedel, martin.raedel@dlr.de
% Copyright (C) 2018-2018 DLR Composite Structures & Adaptive Systems
%                                                               __/|__
%                                                              /_/_/_/  
%                                                                |/ DLR

%
% Package Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{bootstrap_dlrreprt}
  [2018/01/31 v0.1 DLR Bootstrap Report Layout]

%
% Declare class options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Package to process options
\RequirePackage{xkeyval}          % Load key/value package

% frutiger
\define@boolkey{bootstrap_dlrreprt.cls}[bootstrapdlrreprt@]{frutiger}[true]{}
%\ExecuteOptionsX{frutiger}

% layout
\define@boolkey{bootstrap_dlrreprt.cls}[bootstrapdlrreprt@]{elements}[true]{}
\ExecuteOptionsX{elements}

% title
\define@boolkey{bootstrap_dlrreprt.cls}[bootstrapdlrreprt@]{title}[true]{}
\ExecuteOptionsX{title}

% pages
\define@boolkey{bootstrap_dlrreprt.cls}[bootstrapdlrreprt@]{page}[true]{}
\ExecuteOptionsX{page}

% any other options are passed to scrreprt
\DeclareOptionX*{%
  \PassOptionsToClass{\CurrentOption}{scrreprt}%
}

% Process and execute options
\ProcessOptionsX\relax

% Load scrreprt class
\LoadClass{scrreprt}

%
% Load Packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{eso-pic}
\RequirePackage{etoolbox}
\RequirePackage{graphicx}         % Load graphicx package (e.g. \resizebox}
\RequirePackage{hyphenat}
\RequirePackage{ifthen}
\RequirePackage[autooneside,headsepline=1pt]{scrlayer-scrpage}% must be loaded after scrreprt
\RequirePackage{tikz}
\RequirePackage{xcolor}
\RequirePackage{xparse}
\RequirePackage{xspace}           % Load intelligent space package

%
% Packages options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usetikzlibrary{calc,fit,positioning,shapes}

%
% Some extra document info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These commands are created to allow for an easy switch between
% dlrreprt and bootstrap_dlrreprt.

\@ifundefined{subject}{     \def\subject#1{\gdef\@subject{#1} }%
                            \def\@subject{}}{}%
\@ifundefined{keywords}{    \def\keywords#1{\gdef\@keywords{#1} }%
                            \def\@keywords{}}{}
\@ifundefined{restrictions}{\def\restrictions#1{\gdef\@restrictions{#1} }%
                            \def\@restrictions{}}{}
\@ifundefined{author}{      \def\author#1{\gdef\@author{#1} }%
                            \def\@author{}}{}
\@ifundefined{coauthor}{    \def\coauthor#1{\gdef\@coauthor{#1} }%
                            \def\@coauthor{}}{}
\@ifundefined{savedby}{     \def\savedby#1{\gdef\@savedby{#1} }%
                            \def\@savedby{}}{}      
\@ifundefined{approvedby}{  \def\approvedby#1{\gdef\@approvedby{#1} }%
                            \def\@approvedby{}}{}  
\@ifundefined{releasedby}{  \def\releasedby#1{\gdef\@releasedby{#1} }%
                            \def\@releasedby{}}{} 
\@ifundefined{contactname}{ \def\contactname#1{\gdef\@contactname{#1} }%
                            \def\@contactname{}}{}
\@ifundefined{telnumber}{   \def\telnumber#1{\gdef\@telnumber{#1} }%
                            \def\@telnumber{}}{}
\@ifundefined{faxnumber}{   \def\faxnumber#1{\gdef\@faxnumber{#1} }%
                            \def\@faxnumber{}}{}      
\@ifundefined{email}{       \def\email#1{\gdef\@email{#1} }%
                            \def\@email{}}{}
\@ifundefined{reportno}{    \def\reportno#1{\gdef\@reportno{#1} }%
                            \def\@reportno{}}{}

% From dlrinstitutes.sty
\@ifundefined{dlrinstitute}{\def\dlrinstitute#1{\gdef\@dlrinstitute{#1} }%
                            \def\@dlrinstitute{}}{}
\@ifundefined{dlrsite}{     \def\dlrsite#1{\gdef\@dlrsite{#1} }%
                            \def\@dlrsite{}}{}

%
% Dummy commands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From listofversions.sty
\newcommand{\addversion}[4]{}
\newcommand{\currentversion}[1][version]{}
\newcommand{\addversionssvn}[1][]{}

% Define a dummy to avoid errors for \settitle
\newcommand{\settitle}[1]{}

%
% Colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\definecolor{DLRgray}{rgb}{0.27,0.27,0.27}

%
% Logos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eagle for use in Tikz
\DeclareDocumentCommand{\tikzdlradler}{D(){DLReagel} O{(0,0)}}{
  %\fontfamily{pfr}\selectfont % Change to frutiger font
  \upshape\mdseries           % Select right shape and series	
  \node (Eagle1SW) at #2                                 [shape=coordinate]{};
  \node (Eagle1SE) at ($(Eagle1SW)  + (3.1em,  0pt)$)    [shape=coordinate]{};
  \node (Eagle1NW) at ($(Eagle1SW)!.333!45:(Eagle1SE)$)  [shape=coordinate]{};
  \node (Eagle1NE) at ($(Eagle1SE)!.333!225:(Eagle1SW)$) [shape=coordinate]{};

  \node (Eagle2NW) at ($(Eagle1SW)!.333!(Eagle1SE)$) [shape=coordinate]{};
  \node (Eagle2SW) at ($(Eagle2NW)!1!90:(Eagle1SW)$) [shape=coordinate]{};
  \node (Eagle2SE) at ($(Eagle1NE)!.333!(Eagle1NW)$) [shape=coordinate]{};
  \node (Eagle2NE) at ($(Eagle2SE)!1!90:(Eagle1NE)$) [shape=coordinate]{};
  
  \draw [line width = .25em, color=DLRgray] (Eagle1NW) -- (Eagle1SW) -- (Eagle1SE) -- (Eagle1NE) -- cycle
                                                  (Eagle2NW) -- (Eagle2SW) -- (Eagle2SE) -- (Eagle2NE) -- cycle;
                                                  
  \node (EagleDLR1) at ($(Eagle2SW)!1!-90:(Eagle2NW)+(0pt, -.2em)$)        [shape=coordinate]{};
  \node (EagleDLR2) [right=0pt of EagleDLR1, inner sep=0pt, align=left, color=DLRgray]{\bfseries DLR};
  
  % Add a box arround it to refer to the logo using label
  \node(#1) [inner sep=0, fit=(EagleDLR2) (Eagle1SW) (Eagle2NE)]{};
}

% Logo for use in Tikz
\DeclareDocumentCommand{\tikzdlrlogo}{D(){DLReagel} O{(0,0)}}{
  \tikzdlradler(#1)[#2]%
  \node (EagleDLR3)  [right=1em of EagleDLR2, inner sep=0pt, align=left, color=DLRgray]{\bfseries f\"ur Luft- und Raumfahrt};
  \node (EagleDLR3a) [above=1.2em of EagleDLR3.base west, anchor=base west, inner sep=0pt, align=left, color=DLRgray]{\bfseries Deutsches Zentrum};		
  \ifnum\language=0%
    \edef\tikzdlrlogo@check@thelanguge{\string english}%              language 0 is always english (BUG in article.cls)
  \else%
    \edef\tikzdlrlogo@check@thelanguge{\languagename}%                get language
  \fi%
  \edef\tikzdlrlogo@check@german{\string german}%                     define helper to check german
  \edef\tikzdlrlogo@check@ngerman{\string ngerman}%                   define helper to check ngerman
  \ifx\tikzdlrlogo@check@thelanguge\tikzdlrlogo@check@german%         check german
  \else\ifx\tikzdlrlogo@check@thelanguge\tikzdlrlogo@check@ngerman%   check ngerman
  \else%                                                              other language
    \node (EagleDLR3b) [below=1.3em of EagleDLR3.base west, anchor=base west, inner sep=0pt, align=left, color=DLRgray]{\scalebox {0.97}{\mdseries German Aerospace Center}};
  \fi%
  \fi%
}

% Eagle and Logo for use without Tikz
\newcommand{\dlradler}      {\begin{tikzpicture}\tikzdlradler\end{tikzpicture}}
\newcommand{\dlrlogo}       {\parbox{18.15em}{\begin{tikzpicture}\tikzdlrlogo\end{tikzpicture}}}
\newcommand{\minidlradler}  {\raisebox{-.5ex}{\begin{tikzpicture}\tikzminidlradler\end{tikzpicture}}\xspace}
\newcommand{\dlradlerinline}{\raisebox{-0.4em}{\scalebox{0.4}{\dlradler}}\xspace}

% DLR Hook
\newcommand{\dlrhook}{%
  \raisebox{.06em}{%
    \begin{tikzpicture}[rounded corners=0]
      % Define nodes for hook
      \node (corner)  at (0,0)                       [shape=coordinate]{};
      \node (left)    at ($(corner)+(-1em, 0)$)      [shape=coordinate]{};
      \node (below)   at ($(corner)+(-.5em, -.5em)$) [shape=coordinate]{};
      % Define nodes for frame 
      \node (FRoutcorner) at ($(corner)+(.09em, .04em)$)           [shape=coordinate]{};
      \node (FRleftout)   at ($(FRoutcorner)+(-.85em, 0)$)         [shape=coordinate]{};
      \node (FRleftin)    at ($(FRleftout)+(.2em, -.2em)$)         [shape=coordinate]{};
      \node (FRbelowout)  at ($(FRoutcorner)!.6em!45:(FRleftout)$) [shape=coordinate]{};
      \node (FRbelowin)   at ($(FRbelowout)+(-.2em, 0)$)           [shape=coordinate]{};
      % Define frame using clip
      \clip(FRleftin) -- (FRleftout) -- (FRoutcorner) -- (FRbelowout) -- (FRbelowin) -- cycle;
      % Draw clip
      \draw[line width = .1em]  (left) -- (corner) -- (below);
    \end{tikzpicture}
  }
}

%
% Fonts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set to Sans-Serif
% \renewcommand*\familydefault{\sfdefault}

% Set fonts for titles
\setkomafont{chapter}{%           % "Kapitel"
  \fontsize{24.0pt}{29.0pt}%      % 
  \color{black!73}%               % 73% black
  \selectfont
}

\setkomafont{section}{%           % "Überschrift"
  \fontsize{17.28pt}{21.5pt}%     % 
  \color{black!73}%               % 73% black
  \selectfont
}

%
% Layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\setlength{\parindent}{0pt}       % No paragraph indent
\KOMAoption{parskip}{full}        % Skip one line for each paragraph

%
% Provide a language switch according to babel -> not required yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\AtBeginDocument{
%  \@ifpackageloaded{babel}{
%    \let\bootstrap_dlrreprt@iflanguage\iflanguage}{
%    \def\bootstrap_dlrreprt@iflanguage#1#2#3{#3}
%  }
%}

%
% Document info available outside of the maketitle commands from getdocinfo.sty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generically make commands available outside
\def\getdocinfo#1{%                              % Define \getdocinfo{author} to...
          \expandafter\def%                      % ... def
          \csname the#1\expandafter\endcsname%   % ... ... \theauthor
          \expandafter%                          % ... ... =
         {\csname @#1\endcsname}%               % ... ... \@author
  \gpreto{\maketitle}%                           % ... at the beginning of the \maketitle command
         {\expandafter\let%                      % ... let
          \csname the#1\expandafter\endcsname%   % ... \theauthor
          \csname @#1\endcsname}%                % ... = \@author
}

% Make typical commands available
\getdocinfo{title}
\getdocinfo{author}
\getdocinfo{publishers}
\getdocinfo{reportno}
\getdocinfo{subject}
\getdocinfo{keywords}
\getdocinfo{restrictions}
\getdocinfo{coauthor}
\getdocinfo{savedby}
\getdocinfo{approvedby}
\getdocinfo{releasedby}
\getdocinfo{contactname}
\getdocinfo{telnumber}
\getdocinfo{faxnumber}
\getdocinfo{email}
\getdocinfo{date}

%
% Adjust the commands useable by hyperref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\newcommand{\theurl}{\INST[\@currinst]{url}}
\@ifpackageloaded{hyperref}{
  \renewcommand{\theemail}{\href{mailto:\theemail}{\theemail}}
  %\renewcommand{\theurl}{\href{\theurl}{\theurl}}
}

%
% An error command to make sure, the contacts are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\bootstrap_dlrreprt@missing#1{
	\PackageError{bootstrap_dlreprt}
	             {No \@backslashchar#1 defined}
	             {Please define it}}
\author   {\bootstrap_dlrreprt@missing{author}}
\telnumber{\bootstrap_dlrreprt@missing{telnumber}}
\faxnumber{\bootstrap_dlrreprt@missing{faxnumber}}
\email    {\bootstrap_dlrreprt@missing{email}}

%
% Elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifbootstrapdlrreprt@elements      
  
  % Table
  % Format caption on top
  \RequirePackage{scrhack}%        Prevent warnings by float/scr packages
  \RequirePackage{float}%          Restyle floats
  \newcommand\floatc@dlrtable[2]{\centering {#1: \scr@fnt@caption #2}} % Format of caption
  \newcommand\fs@dlrtable{%        Style of dlrtable
    \let\@fs@iftopcapt\iftrue%     ... captions on top
    \let\@fs@capt\floatc@dlrtable% ... caption format as above
    \def\@fs@pre  { }%             ... nothing before
    \def\@fs@mid  { \kern0.5ex }%  ... space under caption
    \def\@fs@post { }%             ... nothing after
  }
  \floatstyle{dlrtable}%           Apply the style dlrtable
  \restylefloat{table}%            ... to the table float

  % create \rowstyle to define a style for a complete row
  \RequirePackage{etoolbox}%                  To allow patching commands
  \newcommand{\rowstyle}[1]{%                 Define new \rowstyle:
    \gdef\currentrowstyle{#1}%                ... set intermediate \currentrowstyle 
    \gdef\resetrowstyle{yes}%                 ... remember to change \\
    #1\ignorespaces}%                         ... run row style once for first cell
  \newcommand{\pretocell}{%                   Define new \pretocell
    \@ifundefined{currentrowstyle}{}{%        ... only if intermediate \currentrowstyle is defined
      \currentrowstyle%                       ... run it
      \@ifundefined{resetrowstyle}{}{%        ... patch \\ commands, if requested:
        \global\let\orig@tnl\tabularnewline%  ...                                      ... save old \tabularnewline
        \global\let\orig@tbs\\%               ...                                      ... and \\ commands
        \gpreto{\tabularnewline}{\global\let\currentrowstyle\@undefined%               ... reset intermediate \currentrowstyle at end line
                                 \global\let\tabularnewline\orig@tnl%                  ... restore original \tabularnewline
                                 \global\let\\\orig@tbs}%                              ... and \\ commands
        \gpreto{\\}             {\global\let\currentrowstyle\@undefined%               ... reset intermediate \currentrowstyle at end line
                                 \global\let\tabularnewline\orig@tnl%                  ... restore original \tabularnewline
                                 \global\let\\\orig@tbs}%                              ... and \\ commands
        \gappto{\tabularnewline}{\addlinespace}%                                       ... add a linespace after line break
        \gappto{\\}             {\addlinespace}%                                       ... add a linespace after line break
        \global\let\resetrowstyle\@undefined% ...                                      ... unremember to patch
  }}}
  \RequirePackage{array}%                     To allow changing column types
  \newcolumntype{l}{>{\pretocell}{l}}%        Add \pretocell to each column type
  \newcolumntype{c}{>{\pretocell}{c}}%        .
  \newcolumntype{r}{>{\pretocell}{r}}%        .
  \newcolumntype{p}{>{\pretocell}{p}}%        .
  \newcolumntype{m}{>{\pretocell}{m}}%        .
  \newcolumntype{b}{>{\pretocell}{b}}%        .
  
  % DLR table head lines: \dlrtablehead
  \RequirePackage{colortbl}%                                       Use colored table backgrounds
  \RequirePackage{booktabs}%                                       Use nice horizontal lines
  \RequirePackage{bigstrut}%                                       Use improved vertical space in cells
  \newcommand{\dlrtablehead}{\specialrule{\heavyrulewidth}{0pt}{0pt}% \dlrtablehead = toprule without distance
                          \rowstyle{\cellcolor{black!15}%                       75% gray in background
                                    \bfseries}%                                 bold font
                          \hspace{0pt}\bigstrut%                                more space in cell
                         }
  % Use DLR hook for first level itemize
  \renewcommand{\labelitemi}{\dlrhook}
\fi

%
% Titlepage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execute title options
\ifbootstrapdlrreprt@title
  \renewcommand{\maketitle}{
  
  \begin{titlepage}
%
      % Make logo appear lower than usual...
      \enlargethispage{3cm}%
%
      % Width of the left/top page margin and the typical text fields
      \newlength{\totalleft}  \setlength  {\totalleft}  {\oddsidemargin}%
                              \addtolength{\totalleft}  {1in}%
                              \addtolength{\totalleft}  {\hoffset}%
      \newlength{\totaltop}   \setlength  {\totaltop}   {\topmargin}%
                              \addtolength{\totaltop}   {1in}%     
                              \addtolength{\totaltop}   {\headheight}%     
                              \addtolength{\totaltop}   {\headsep}%     
                              \addtolength{\totaltop}   {\topskip}%     
                              \addtolength{\totaltop}   {\voffset}%     
                              \addtolength{\totaltop}   {14.8pt}%
      \newlength{\fieldindent}\setlength  {\fieldindent}{0.366\paperwidth}%
      \newlength{\fieldmargin}\setlength  {\fieldmargin}{1cm}%
      \newlength{\fieldwidth} \setlength  {\fieldwidth} {\paperwidth}%     
                              \addtolength{\fieldwidth} {-\fieldindent}%     
                              \addtolength{\fieldwidth} {-\fieldmargin}%
%
      % Create background picture  
      \AddToShipoutPicture*{%                                              % Add the background to shipout
        \resizebox{\paperwidth}{\paperheight}{%                            % Resize to page size
          \begin{tikzpicture}%                                             % Draw a tikzpicture!
            \useasboundingbox (   0mm,   0mm) rectangle (210mm,297mm);%    % The bounding box of the picture
            \clip             ( 184mm, 618mm) circle    (405mm);%          % Clip "horizontal" edge
            \clip             ( 447mm, 369mm) circle    (405mm);%          % Clip "vertical" edge
            \fill[black!73]   (-001mm,-001mm) rectangle (211mm,298mm);%    % Fill the remaining area with gray/silver
          \end{tikzpicture}                                                %. end tikzpicture
        }%                                                                 %. end resizebox
      }%                                                                   %. end background
%
      % Font settings
      \sffamily%                       % ... and activate it
      \sloppy%                         % Don't break words
      \hfuzz=5cm%                      % Don't warn for bad boxes
%
      % Move cursor to topleft corner
      \vspace*{-\totaltop}\Large%
      \begin{addmargin}[-\totalleft]{0cm}%
        \begin{tikzpicture}
          \node (hook) at (0,0) {};
          \node (subject)
              at (\fieldindent,-0.12\paperheight) 
                 [anchor         = north west,
                  text width     = \fieldwidth, 
                  minimum height = 0, 
                  inner sep      = 0,
                  align          = left]
                 {%
                   \color{white} {\Huge \bfseries \@subject} \\[.5ex]% The subject
                   %{\LARGE \bfseries \@reportno}% The report nr. -> I choose not to print this, as it problably not of interest for the "outside" world. MR
                 };
          \node (title)
              at (\fieldindent,-0.37\paperheight) 
                 [anchor         = north west,
                  text width     = \fieldwidth, 
                  minimum height = 0, 
                  inner sep      = 0,
                  align          = left]
                 {\color{black!73} \huge \bfseries \nohyphens{\@title}};
          \node  (info)
                 [below=2em of title.south west, 
                  anchor=north west,
                  text width = \fieldwidth, 
                  minimum height=0, 
                  inner sep =0,
                  align = left]
                 {\color{black!73} \@subtitle\\[1.3em]
                                   \@author\\[1.3em]
                                   \@publishers};
        \end{tikzpicture}
        \vfill\large%
        \begin{tikzpicture}
          \node (hook) at (0,0) {};
          \tikzdlrlogo[(\fieldindent,-2cm)]
        \end{tikzpicture}
      \end{addmargin}
    \end{titlepage}
    \cleardoublepage
  }% End renewcommand
  
\fi

%
% Header/Footer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifbootstrapdlrreprt@page
  % Fonts for headers/footers
  \setkomafont{pagefoot}{\sffamily\itshape}
  
  % Select headers/footers style
  \pagestyle{scrheadings} 
  
  %Clear old  headers/footers
  \clearscrheadfoot
  
  % Define chapter and section title for header
  \automark[section]{chapter}       % Kapitel (und Section) in Kopfzeile
  
  % Create new header
  \makeatletter
    \if@twoside
      \rehead[]{\headmark}
      \lohead[]{\headmark}
    \else
      \ihead[]{\headmark}
    \fi 
  \makeatother
  \ohead[]{\pagemark}
  
  % Create new footer
  \makeatletter
  \newcommand{\copyrighttext}{
    {%
      \dlradlerinline%
      %DLR%
    } -- \@title%
  }
  \makeatother

  \ofoot[\pagemark]{}
  \cfoot[\copyrighttext]{\copyrighttext}
\fi

% eof
