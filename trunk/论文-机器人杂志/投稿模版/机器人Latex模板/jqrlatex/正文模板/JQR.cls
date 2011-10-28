%%
%% This is file `JQR.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% JQR.dtx  (with options: `cls,art')
%% 
%% This is a generated file.
%% 
%% Copyright 2005-2006
%% CTEX.ORG and any individual authors listed elsewhere in this file.
%% 
%% This file was generated from file(s) of the JQR package project.
%% -----------------------------------------------------------------
%% 
%% It may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License, either version 1.2
%% of this license or (at your option) any later version.
%% The laJQR version of this license is in
%%    http://www.latex-project.org/lppl.txt
%% and version 1.2 or later is part of all distributions of LaTeX
%% version 1999/12/01 or later.
%% 
%% This file may only be distributed together with a copy of original
%% source files. You may however distribute original source files
%% without such generated files. Copying of this file is authorized
%% only if either:
%% (1) you make absolutely no changes to your copy, including name; OR
%% (2) if you do make changes, you first rename it to some other name.
%% 
%% To produce the documentation run the original source files ending
%% with `.dtx' through LaTeX.
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{JQR}
  [2005/11/28 v0.1a JQR
   article document class]
\newif\ifJQR@dvips \JQR@dvipsfalse
\DeclareOption{dvips}{\JQR@dvipstrue}
\DeclareOption{dvipdfm}{\JQR@dvipsfalse}
\newif\ifJQR@english \JQR@englishfalse
\DeclareOption{english}{\JQR@englishtrue}
\DeclareOption{chinese}{\JQR@englishfalse}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{ctexart}}
\ProcessOptions
\LoadClass[a4paper,twoside,fancyhdr]{ctexart}[2004/08/14]
%% math packages
\RequirePackage{amsmath,amsthm,amsfonts,amssymb,multicol,times,mathptmx,bm,array,multirow,wasysym,slashbox}
%% graphics packages
\RequirePackage{graphicx}

\DeclareSymbolFont{myletters}{OML}{cmm}{m}{it}
\DeclareMathSymbol{\myepsilon}{\mathord}{myletters}{"0F}

\AtEndOfPackage{\makeatletter\input{JQR.cfg}\makeatother}
%%页面设置%%
\headheight 7 true mm 
\topmargin 0pt 
\oddsidemargin 0pt
\evensidemargin 0pt 
\textheight 245 true mm 
\textwidth 173 true mm
\marginparwidth 8 true mm 
\columnsep 8 true mm 
\footskip 250pt
\headsep 3mm 
\hoffset -8mm 
\voffset -5mm
\parskip 0pt
\parindent=2\ccwd
\abovedisplayskip=3pt plus 0pt minus 1pt
\belowdisplayskip=3pt plus 0pt minus 1pt
\jot=1pt
\gdef\arraystretch{1.16}\selectfont
\gdef\baselinestretch{1.233}\selectfont
\bibliographystyle{plain}
%%
%% the chinese or english environments
%%
\def\chinesematter{%
  \def\baselinestretch{1.2}\selectfont
  \def\JQR@abstract{\qquad摘\quad 要：}
  \def\JQR@keywords{\qquad关键词：}
}
\def\englishmatter{%
  \def\baselinestretch{1.2}\selectfont
  \def\JQR@abstract{\quad Abstract: }
  \def\JQR@keywords{\quad Keywords: }
}
%%
%% the section formats
%%
\def\CTEX@section@format{\bfseries\zihao{-4}}
\def\CTEX@section@indent{\z@}
\def\CTEX@section@beforeskip{-1.2ex \@plus -.5ex \@minus -1ex}
\def\CTEX@section@afterskip{.5ex \@plus .1ex \@minus .1ex}

\def\CTEX@subsection@format{\bfseries\zihao{5}}
\def\CTEX@subsection@indent{\z@}
\def\CTEX@subsection@beforeskip{\z@}
\def\CTEX@subsection@afterskip{.1ex \@minus .1ex}

\def\CTEX@subsubsection@format{\bfseries\zihao{5}}
\def\CTEX@subsubsection@indent{\z@}
\def\CTEX@subsubsection@beforeskip{\z@}
\def\CTEX@subsubsection@afterskip{.1ex \@minus .1ex}

\gdef\@author{}
\renewcommand*{\author}[2][\@empty]{%
  \ifx \@author\@empty
    \gdef\@author{#2\kern 0.5pt \textsuperscript{#1}}
  \else
    \expandafter\gdef\expandafter\@author\expandafter{%
      \@author，{#2\kern 0.5pt \textsuperscript{#1}}}
  \fi}
\gdef\@address{}
\newcommand*{\address}[2][\@empty]{%
  \ifx \@address\@empty
    \gdef\@address{#1#2}
  \else
    \expandafter\gdef\expandafter\@address\expandafter{%
      \@address\tabularnewline{#1#2}}
  \fi}
\renewcommand{\thefootnote}{}
\gdef\thefootnote{}\nofiles
\renewcommand\maketitle{\par
  \begingroup
      \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\@title\@empty
  \global\let\@author\@empty
  \global\let\@address\@empty
  \global\let\@thanks\@empty
}
\newenvironment{chinesetitle}
  {\chinesematter
   \def\@maketitle{%
  \begingroup\centering%
   {\bfseries\zihao{3} \@title \par}%
    \vskip 3mm%
    {\fangsong\large\zihao{4}
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    {\normalfont\normalsize\zihao{6}
        （\@address）
      \par\vskip 18pt}%
  \endgroup}}
  {\vspace*{7mm}}
\newenvironment{englishtitle}
  {\englishmatter
   \def\@maketitle{%
  \begingroup\centering%
   {\bfseries\zihao{-4} \@title \par}%
    \vskip 3mm%
    {\fangsong\large\zihao{5}
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    {\normalfont\normalsize\zihao{6}
        (\@address)
      \par\vskip 18pt}%
  \endgroup}}
  {\vspace*{6mm}}
\AtBeginDocument{\null}
\renewenvironment{abstract}
  {\begin{minipage}{158mm}\zihao{-5}%
   \textbf{\JQR@abstract}\ignorespaces}
   {\end{minipage}}
  \newcommand*\keywords[1]{%
  \textbf{\JQR@keywords}\begin{minipage}[t]{139mm}#1\end{minipage}\par}
  \newcommand{\category}[2]{%
  \zihao{-5}\qquad\textbf{中图分类号：}#1\kern 20mm\textbf{文献标识码：}#2}
  \newcommand*\doi[1]{\noindent\zihao{-5}
  \textbf{DOI}：10.3724/SP.J.1218.2010.00#1\par\vskip 8mm\zihao{5}}
  \renewcommand{\footnoterule}{\vspace*{2mm}\vbox to 0pt {\hbox to 5cm{\hrulefill}\vss}\vspace{3mm}}
  \newcommand{\dates}[2]{%
   \footnotetext{\zihao{6}通讯作者：#1\qquad 收稿\!／\!\!录用\!／\!\!修回：#2}\par\footnotetext{}}
  \newcommand{\supported}[1]{%
  \footnotetext{\zihao{6}\hangafter 1\hangindent 6.7em基金项目：#1\kern -2pt．}}
  \def\Intro #1（#2{%
   \hangafter 1\hangindent 4.5em\noindent\quad\zihao{-5}{\kaishu #1}\rm\songti（#2}
   \newenvironment{AuthorIntroduction}
  {\noindent\bfseries\zihao{5}作者简介：\zihao{-5}\normalfont\par\vspace*{1mm}}
  {}
\setcounter{secnumdepth}{3}
\newtheoremstyle{break}% name
  {}%      Space above, empty = `usual value'
  {}%      Space below
  {\itshape}% Body font
  {}%         Indent amount (empty = no indent, \parindent = para indent)
  {\bfseries}% Thm head font
  {.}%        Punctuation after thm head
  {\newline}% Space after thm head: \newline = linebreak
  {}%         Thm head spec
\renewenvironment{thebibliography}[1]
     {\fontsize{8.5}{9.8}\normalfont
      \list{\@biblabel{\@arabic\c@enumiv}\hfill}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \itemsep -2\p@ \@plus 1\p@ \@minus 1\p@
      \parsep 0ex
      \topsep 0ex
      \partopsep 0ex
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newenvironment{References}
 {\vskip 1ex \centering \zihao{5}\bfseries{参考文献 （References）}\vskip .5ex%
 \begin{thebibliography}{99}}
 {\end{thebibliography}}
\RequirePackage{ifthen}
\makeatletter
\newcommand{\allmodesymb}[2]{\relax\ifmmode{\mathchoice
{\mbox{\fontsize{\tf@size}{\tf@size}#1{#2}}}
{\mbox{\fontsize{\tf@size}{\tf@size}#1{#2}}}
{\mbox{\fontsize{\sf@size}{\sf@size}#1{#2}}}
{\mbox{\fontsize{\ssf@size}{\ssf@size}#1{#2}}}}
\else
\mbox{#1{#2}}\fi}
\makeatother
\newcommand{\greeksym}[1]{\usefont{U}{psy}{m}{n}#1}
\newcommand{\ualpha}{\allmodesymb{\greeksym}{a}}
\newcommand{\ubeta}{\allmodesymb{\greeksym}{b}}
\newcommand{\uchi}{\allmodesymb{\greeksym}{c}}
\newcommand{\udelta}{\allmodesymb{\greeksym}{d}}
\newcommand{\ugamma}{\allmodesymb{\greeksym}{g}}
\newcommand{\umu}{\allmodesymb{\greeksym}{m}}
\newcommand{\unu}{\allmodesymb{\greeksym}{n}}
\newcommand{\upi}{\allmodesymb{\greeksym}{p}}
\newcommand{\utau}{\allmodesymb{\greeksym}{t}}
\newcommand{\utheta}{\allmodesymb{\greeksym}{q}}
\endinput
%%
%% End of file `JQR.cls'.
