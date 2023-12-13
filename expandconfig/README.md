# expandconfig
Small program which reads text files where there can be dynamic text
parts with several alternatives and returns text files with all
possible combinations.

## Backstory
I have to make numerical experiments using a program which reads config
files.
Like in all numerical experiments, some parameters in the config file
have to varied.
Instead of manually creating hundreds of config files I decided to
write a small program which does this for me.

Fortunately, this program works not just on config files but on all
text files (at least I guess so).

## Examples
### Basics
Consider the config file `basic.cfg`:

    this = [[ 1   || 2    ]]
    that = [[ 0.5 || 0.75 ]]

A call `expandconfig basic.cfg` will return four different files
`basic.cfg.0000001.exp` up to `basic.cfg.0000004.exp` with the contents

    this =  1   
    that =  0.5 

and

    this =  2    
    that =  0.5 

and

    this =  1   
    that =  0.75  

and finally

    this =  2    
    that =  0.75 

### Tagging
It is possible to give tags to alternative text parts.
Parts with the same tag will not be varied independent of each other,
but rather together.
Consider the config file `tags.cfg`:

    this = [[((foo)) 1   || 2    ]]
    that = [[((foo)) 0.5 || 0.75 ]]

A call `expandconfig tags.cfg` will return only two different files
`tags.cfg.0000001.exp` and `tags.cfg.0000002.exp` with the contents

    this =  1   
    that =  0.5 

and

    this =  2    
    that =  0.75 

Tags will never be printed.
Of course, dynamic parts without tags can be mixed with dynamic
parts with tags.

## "Language"

There are no special strings, except for `[[`, `]]`, `||`, `((`
and `))`.
A file that can be parsed with `expanconfig` can be described in
the following way:

    "static part"  [[ (( "tag name" )) "dynamic part" ]]  "static part"
                      \______________/
                             |
                          optional
                   \__________________________________________________/
                                             |
                           can be repeated as often as you want*

where "dynamic part" consists of

    "alternative" || "alternative"
                  \______________/
                         |
          can be repeated as often as you want*

Notice that the spaces are only included due to readability.

The "static part", "tag name" and "alternative" can be almost any
text you want except these few special rules:

 - "static part" must not contain literal `[[`
 - "tag name" must not contain literal `))`
 - "alternative" must not contain literal `||` and `]]`

While (except for these rules) every part can contain even the otherwise
special strings `[[`, `]]`, `||`, `((` and `))`, I stronly discourage
you to use them, because it is rather error-prone.
Otherwise `expandconfig` may not agree with you on what it should do.

The asterisk on the comment that something can be repeated as often as
you want is due to the fact that expandconfig will not create any
files if your input file has more than `9999999 = 1e7 - 1`
combinatorial possibilities.

## Building
Building was tested only on Linux. You need `make`, `gcc`. For
building the project type `make`, for running the test type `make test`,
which will most likely only tell you that the test was passed.

## Bugs?
Everything works as intended at the moment.

## Future improvements?
 - Further improving the quality of the code.
 - Adding (better) error handling.
 - Additional features.

