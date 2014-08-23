// http://mingw-users.1079350.n2.nabble.com/Query-regarding-offered-alternative-to-asprintf-tp6329481p7578880.html

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

int asprintf( char **, char *, ... );
int vasprintf( char **, char *, va_list );
FILE *ftrash = NULL;

int vasprintf( char **sptr, char *fmt, va_list argv )
{
        if(!ftrash) ftrash = fopen("nul", "wb"); /* on windows, /dev/null = nul */
        if(!ftrash) fprintf(stderr, "this shouldn't happen\n");

/* NOTE: remember to fclose(ftrash) after the last call to vasprintf() has been made */

/****
old code:
        int wanted = vsnprintf( *sptr = NULL, 0, fmt, argv );

new code: */
        int wanted = vfprintf( ftrash, fmt, argv );
/****/

        if( (wanted < 0) || ((*sptr = malloc( 1 + wanted )) == NULL) )
        return -1;


/****
old code:
        return vsprintf( *sptr, fmt, argv );

new code: */
        (*sptr)[wanted] = '\0';
        return vsnprintf( *sptr, wanted, fmt, argv );
/****/
}

int asprintf( char **sptr, char *fmt, ... )
{
        int retval;
        va_list argv;
        va_start( argv, fmt );
        retval = vasprintf( sptr, fmt, argv );

        va_end( argv );
        return retval;
}
