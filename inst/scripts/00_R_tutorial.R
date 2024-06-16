#' @title Introduction to R Studio
#'
#' @author Dario Fanucchi, 2012 @author Juliet R. C. Pearson, 2014 @author Carl
#' A. B. Pearson, 2024
#'
#' International Clinics on Infectious Disease Dynamics and Data Program
#'
#' @license Some Rights Reserved, CC BY-NC 4.0
#' (https://creativecommons.org/licenses/by-nc/4.0/) Last updated by Carl A. B.
#' Pearson, June 2024
#'
#' @description RStudio is an Integrated Development Environment (IDE) for the R
#' programming language. RStudio combines for editing source files, running
#' lines or sections of code in an R console, typing into the console directly,
#' working with the file system, working in the terminal, and many more useful
#' tools for coding work, all combined in a single application. RStudio also
#' offers some R-specific behavior, such as keeping track of the variables
#' currently in the workspace and showing the command history.
#'
#' This tutorial introduces the RStudio IDE, which will help you complete the
#' other R tutorials more effective. Once you've worked through the tutorial,
#' you should be able to:
#'
#'  - understand the main tabs in RStudio (Editor, Console, Environment,
#'    History, Help, Plots, Packages)
#'  - run code from the console, either directly or quickly from the editor
#'  - browse the installed packages and load libraries
#'  - make use of the Help system
#'  - use shortcut keys
#'
#' If you are not reading this file from within RStudio, you can open it in
#' RStudio in one of several ways:
#'
#'  - Open the file using "File" -> "Open File" in the RStudio editor
#'  - At the terminal/cmd write the following: rstudio introRstudio.R (or, on a
#'  Mac: open -a RStudio introRstudio.R)
#'  - In RStudio in the files window (bottom right of the screen or Ctrl+5 to
#'  activate), browse to the location of the file by clicking the "..." to the
#'  top-right of the window, and then open the directory. Click on the file in
#'  the files window and it will load in the editor.
#'
#' Note that - throughout the tutorial - Mac users should use the command (open
#' apple) key if the keystroke sequence calls for the control key but does not
#' appear to function as expected.
#'
#' @section Using RStudio to run line-by-line from scripts
#'
#' In RStudio you can execute commands line by line or in the console. Copy the
#' following command (without the leading "#'") into the console and press
#' enter:
#'
#' @example plot(c(1,2,3,4,5), c(1,2,3,4,5))
#'
#' A plot should appear in the plots window as a result of the command. Instead
#' of typing directly in the console you can type the command in a .R script
#' like this one and run it in the console. To do so: place your cursor on the
#' line below and press Ctrl+Enter

plot(c(1, 2, 3, 4, 5), c(1, 2, 3, 4, 5), "l")

#' Notice how a plot appears in the plot window corresponding to this command.
#' Don't worry about exactly how the plot command works yet. Anything that you
#' run from a script using Ctrl+Enter gets typed directly in the R console. You
#' can equivalently copy everything out in the console and it will work fine,
#' but it is generally easier to work in the editor and use Ctrl+Enter to
#' execute commands in the console.
#'
#' @section Variables, the Environment, and basic arithmetic
#'
#' In R you can assign a value to a variable by using the "<-" operator.
#' Run the following line of code by pressing Ctrl+Enter

x <- 10

#' Make sure the Environment window is visible by typing Ctrl+8. Notice that
#' the Environment window now contains a variable `x`, with a value of 10. To
#' remove x from the workspace use the `remove` function:

remove(x)

#' Notice how x no longer shows up in the Environment window. Typing `x` now
#' produces an error. Try the command below in the console:

x  # will produce an "Error: object 'x' not found"!

#' You may also see the `rm()` function used in the same way; it is identical to
#' using `remove()`.
#'
#' R can be used to perform arithmetic in a similar fashion to a calculator.
#' For instance, you can type this:

(1242 - 241.1) * 32.21

#' When you press Ctrl+Enter, the cursor automatically moves to the next line in
#' the editor. This way, you can run scripts one line at a time by repeatedly
#' pressing Ctrl+Enter.
#'
#' The following code does some basic arithmetic.  Run it line by line:

x <- 2              # pick a number between 1 and 10
y <- x * 9          # multiply it by 9
d1 <- floor(y / 10) # get the first digit (floor(x) gives the integer part of x)
d2 <- y %% 10       # get the second digit (a %% b gives the remainder)
d1 + d2 - 4         # add the digits and subtract 4: should always be 5!

#' The code above performed a classic arithmetic "teaser". Starting with any
#' number `x`. Multiply it by 9, add the digits of the result, subtract 4 and
#' the answer will always be 5. Try changing the initial value of x to some
#' other number between 1 and 10 and see if it works. Also notice how it is
#' possible to add comments next to code with a "#" after the code to explain
#' what it does.
#'
#' As you ran the code lines above, the variables `d1`, `d`, `x` and `y` were
#' added to the workspace with their values. You can see in the Environment that
#' `d1` and `d2` are indeed the first and second digit of `y` respectively.
#'
#' Since we're done with these variables, we'll keep our workspace clean by
#' removing them:

rm(x, y, d1, d2) # remember: `rm()` is equivalent to `remove()` used above

#' We can do a *full* wipe of all variables in the Environment using the broom
#' icon in the Environment window, but generally you should prefer to remove
#' objects in a targetted way.
#'
#' @section Running sections of a script
#'
#' Before we begin this next section, you may want to clear your console so that
#' things are less cluttered. Simply type Ctrl+L, and the console will be
#' cleared. If you wish to view your previous commands, you can still place the
#' cursor in the console and press the up arrow to see them. You can press
#' Ctrl+Up-arrow (with your cursor in the console) to see your command history
#' in a pop-up window. If you don't want to reach for a mouse, you can switch
#' back to the editor from the console with Ctrl+1
#'
#' Sometimes you don't want to run your code line-by-line but rather in a batch.
#' In RStudio you can highlight code in the Editor you want to run, press
#' Ctrl+Enter like when running a single line, and now everything highlighted
#' will be run. Let's try it!
#'
#' The code below computes the prime numbers in the range 1-100 and writes them
#' in the console. There are many control structures and commands in this code,
#' and it might feel excessive to run it line-by line. Don't worry about how it
#' works for the moment (although there are line-by-line comments for the very
#' enthusiastic reader). In order to run it all in one shot, highlight the whole
#' section and press Ctrl+Enter.

x <- logical(100)  # make a logical vector of 100 elements, all initial `FALSE`
x[-1] <- TRUE      # Set all but 1 to `TRUE`
prim <- list()     # Setup an Empty list for found primes
for (i in 2:100) { # go through the list from 2 to 100 (we know 1 isn't prime)
  if (x[i] == TRUE) {                  # if i is known to be prime...
    prim[[length(prim) + 1]] <- i      # ... add it to the list of known primes
    if (i <= 50) {                     # ... if i has multiples in [1,100] ...
      for (j in i * (2:(100 / i))) {
        x[j] <- FALSE   # ... and set all multiples of i to be NOT prime
      }
    }
  }
}
prim <- as.integer(prim) # convert list into an integer vector
print(prim) # print the list of discovered primes to the console

#' Notice that all the code gets "typed" into the console when you press
#' Ctrl+Enter. The list of primes below 100 is produced. If you run it one line
#' at a time, it still works, of course - but also note the whole block
#' expressions (that is, part starting from `for (i in `...) are run even though
#' they are multiple lines.
#'
#' @section Vectors in R
#'
#' Look at your Environment (Ctrl+3 if it isn't showing). There are two
#' variables `i` and `j` that are single valued variables similar to those we've
#' seen before. `i` has a value of `100L`. The "L" just refers to the data type
#' associated with `i`: in this case it is a "long integer". There are also two
#' other variables in your workspace that are not as familiar: `x` and `prim`.
#' These are vectors. Next to prim is a summary of the structure of the data it
#' contains: `int [1:25]`. This says that prim is an int(eger) vector with 25
#' values.
#'
#' We can make a new vector by using the `c` (short for "concatenate") function:

myvect <- c(1.1, 2.2, 3.3, 4.4, 5.5, 4.4, 3.3, 2.2, 1.1)

#' Notice that `myvect` is summarized as `num [1:9]`. It contains 9 elements,
#' each of which is a numeric data type.
#'
#' We now have quite a few variables hanging around in the workspace. Clean up
#' the Environment by clicking the "Clear All" button, which looks like a broom,
#' or typing the following code:

rm(list = ls()) # `ls()` gives us the names for all variables in the environment

#' @section Some RStudio-specific tips and tricks
#'
#' @subsection Using code completion and Help
#'
#' Code Completion is an extremely useful feature included in RStudio. If you
#' begin typing the name of some function, you can get RStudio to complete it
#' for you (or provide a list of possible completions) by simply hitting the TAB
#' key. If several completions are possible, they are displayed in a list next
#' to the cursor. If there are options, hitting TAB a second time automatically
#' complete using the first option available; you may instead arrow up / down to
#' select desired completion, then hit TAB to use it.
#'
#' Try completing some of the following phrases by hitting TAB. Note that this
#' works both in the console and the Editor.

plo
do.
len
strs
med

#' Note that if there is only one option available, hitting TAB will
#' automatically complete the command (as with `strs` above).
#'
#' Notice how the autocompletion also gives you a bit more information about the
#' functions. If you need still further information, you can make use of the
#' Help option: Press F1 while highlighting a valid R function or command. For
#' example, try TAB (to get completion list) + F1 of the last statement above:
#' you see the Help page for `median` open.
#'
#' Similarly, press F1 while the cursor is on the line of code below (fn+F1 on a
#' Mac):

median

#' Notice that the Help window (by default in the bottom right corner) loads the
#' help for "Median Value". This feature is very useful indeed, and it can
#' increase your productivity in R substantially to be able to obtain help very
#' quickly on functions you wish to use in this fashion.
#'
#' You also use the following syntax in the console to open associated Help:

?median

#' Lastly, using TAB in an already completed function will show you the
#' arguments for that function and any documentation about them. For example,
#' copy this into the console, put the cursor between the `()`s then hit TAB:
#' you should see the `x=` and `na.rm=` arguments, as well as their usage. Note:
#' this does NOT work in the Editor.

median()

#' @subsection Console History and History tab
#'
#' In addition to getting command history with the arrows in the console, it is
#' also possible to obtain previous statements typed by opening the History tab.
#' The History window can be activated by typing Ctrl+4, or by clicking on the
#' History tab. Inside the History window, you can scroll through all code
#' previously run in the console.
#'
#' Double clicking an entry in the History window (or clicking Enter on it) will
#' cause it to copy to the console, where it can be evaluated. Alternatively,
#' you can click the "To Source" button (or type Shift+Enter) to send the
#' highlighted line(s) to the current location in the source editor.
#'
#' As an example, run the following lines of code:

x <- 1
y <- 1

#' and now run the following block of code

y <- y * (x + 1)
x <- x * (x %% 2 + 1) + 1

#' Place the cursor in the source below this comment, then look in the History
#' window for the above two lines, highlight them and click Shift+Enter. Observe
#' how the lines get copied below. You can type this several times to generate
#' several copies.
#'
#' You can then use Ctrl+1 to get back to the source editor and run the code.
#' Alternatively you could send the code directly to the console from the
#' History. Try both.

# ^ ^ ^ Place cursor above this line when inserting from history ^ ^ ^

#' @section Including Packages
#'
#' To load an installed package, you can go to the Packages tab (Ctrl+7), and
#' select the package that you want to use. For instance, you can include the
#' `Matrix` package in this fashion. You can also do this in code as follows:

require("Matrix")

#' once you've included a package, you can make use of functions from that
#' package. You can click the name of a package in your Packages window to get
#' help for that specific package. Here are some examples using the `Matrix`
#' package

mat <- Matrix(0, 2, 2, sparse = FALSE)
mat[1, 1] <- 2; mat[1, 2] <- 3
mat[2, 1] <- 3; mat[2, 2] <- 4
mat
det(mat)

#' @section Running a whole script
#'
#' There are a few more options for running code in RStudio other than typing
#' Ctrl+Enter. They are provided here with their shortcuts for completeness.
#'   1. Run the ENTIRE script:                              Ctrl+Shift+O
#'   2. Run the script FROM THE BEGINNING UP TO THIS LINE:  Ctrl+Shift+B
#'   3. Run the script FROM THIS LINE TO THE END:           Ctrl+Shift+E
#' These can all be useful in different contexts. You can try them out by moving
#' to various places in this script and typing these key combinations. (Note,
#' however, that there is an error produced on line 96, so you will not be able
#' to run the entire script in this case.)
#'
#' For additional keyboard shortcuts, see Tools > Keyboard Shortcuts Help.
#'
#' @section Appendix. List of useful RStudio Shortcut keys
#' [For a COMPLETE list of keyboard shortcuts type Alt+Shift+K (on Linux
#' or Windows) or Option+Shift+K (on a Mac). Note that some of the shortcuts
#' below work differently on Macs.]
#'
#' Ctrl+1          (Move cursor to Editor)
#' Ctrl+2          (Move cursor to Console)
#' Ctrl+3          (Show Environment window)
#' Ctrl+Enter      (Editor: run a line, or a highlighted section of code)
#' Ctrl+L          (clear the console)
#' Ctrl+Up-arrow   (Console: see list of previous commands; Editor: go to top of
#' file; use Cmd+Up-arrow on a Mac)
#' Ctrl+Down-arrow (Editor: Go to bottom of file; use Cmd+Down-arrow on a Mac)
#' Ctrl+Shift+R    (Editor: Run the entire script)
#' Ctrl+Shift+B    (Editor: Run the script from beginning to this point)
#' Ctrl+Shift+E    (Editor: Run the script from this point to end)
#' F1              (Editor: Get help on whatever is under the cursor)
#' F2              (Editor: Show source code for whatever is under the cursor -
#' if it exists)
#' Ctrl+/          (Editor: Toggle Comment on selected region)
