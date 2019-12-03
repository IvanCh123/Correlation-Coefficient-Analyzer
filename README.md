# Correlation Coefficient Analyzer

## General Description
The Correlation Coefficient Analyzer is a command-line program which provides a bioinformatic-oriented tool to determinate the correlation between different types of cancer, using as criteria the similarity in their genetic content. It's required that the user enters a [.csv](https://en.wikipedia.org/wiki/Comma-separated_values) or [.tsv](https://en.wikipedia.org/wiki/Tab-separated_values) file, whose cells contain normalized floating-point values that are interpreted as the expression of each gene (rows) in a specific type of cancer (columns). The expected input could be observed as follows:

<center>
    <img src="https://i.imgur.com/cc6WATo.png" width="500" height="300">
</center>

This program will generate a reduced version of the initial input file, using the [Pearson Correlation Coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) as a mechanism for discarding data. The user indicates whether reduction is made based on the correlation or anticorrelation, and their ranges, if not specified, these will be set by default as [0.75 , 1]  and  [-1 , -0.75] respectively. The program must validate these ranges for each of the cases and, if necessary, it will generate an error message.


The typical way of operation is the following: the correlation coefficient between a specific type of cancer and each one of the rest is calculated, if it doesn't match the criteria established by the user, for none of the cases, the corresponding column will be discarded. This process is repeated for each of the columns, until the maximum reduction of the table is obtained.

Finally, wildcards can be used as a simplified way to request the program to operate on a specific set. For example, the user can indicate that all types of cancer that contain the word 'carcinoma' are compared with the rest.

## User Manual

As mentioned above, the program will use the console for all its commands, the generic command to use will be ``corr fileName.(csv/tsv) [opts] [Regex]``, where ``[opts]`` is the different options of flags that the program has and ``[regex]`` is an optional space that will only be used if the user wants to make comparisons in a specific set.

#### Command line flags (opts) and their use:

> ``--help`` Displays all the flags and their explanation.
>
> ``-t`` Transposes the matrix.
> 
> ``-cc x:y`` The program will only take into account the correlation between different types of cancer. 
>
> ``-ac x:y`` The program will only take into account the anticorrelation between different types of cancer. 
> 
> ``-o outputFile.csv/tsv`` Exits the program and generates a ``.csv`` or ``.tsv`` file with the summarized table.


* The ``x:y`` in the flags is the range to use. When it comes to correlation, these two digits can only be positive numbers, if the user doesn't indicates them, the default to use will be [0.75, 1]. On the other hand, when it comes to the anticorrelation these two numbers must  be negative, if not, the default range will be [-1, -0.75].

## Output

Once the user types the ``-o`` flag, the program will end and will generate, whether a ``.csv`` file or a ``.tsv`` file (the user chooses this), with the initial table but summarized. The way the program summarizes is the following: if the user used the ``-cc`` flag, the cancers left  in the final file will only be the ones that are correlationated between them, if the user used the ``-ac`` is the opposite. Hence, if a cancer dissapears of the table, is because that cancer doesn't has any correlation/anti with any of the others cancers in the file.
    
### Resources

* [C Programming Language](http://www.cplusplus.com/reference/clibrary/)

* [Command-line Interface](https://en.wikipedia.org/wiki/Command-line_interface)


### About this project

This is an assigment of the Parallel and Concurrent Programing course at the University of Costa Rica, under the guidance of Jeisson Hidalgo Céspedes, and is developed by:

* Gloriana Mora Villalta. Email: gloriana.moravillalta@ucr.ac.cr.

* Iván Chavarría Vega. Email: ivan.chavarriavega@ucr.ac.cr.

This project is protected under the [GNU](https://www.gnu.org/licenses/gpl-3.0.en.html) General Public License.

