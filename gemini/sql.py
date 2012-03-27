# simpleSQL.py
#
# simple demo of using the parsing library to do simple-minded SQL parsing
# could be extended to include where clauses etc.
#
# Copyright (c) 2003, Paul McGuire
# Modified by Aaron Quinlan, 2012
#
from pyparsing import Literal, CaselessLiteral, Word, Upcase, delimitedList, Optional, \
    Combine, Group, alphas, nums, alphanums, ParseException, Forward, oneOf, quotedString, \
    ZeroOrMore, restOfLine, Keyword

# define SQL tokens
selectStmt = Forward()
selectToken = Keyword("select", caseless=True)
fromToken = Keyword("from", caseless=True)

# ARQ 2012-Feb-10: allow struct-like column names, e.,g. gt_types.sample1 (add + ".$")
ident = Word( alphas, alphanums + "_$" + ".$" ).setName("identifier")

columnName = Upcase( delimitedList( ident, ".", combine=True ) )
columnNameList = Group( delimitedList( columnName ) )
tableName = Upcase( delimitedList( ident, ".", combine=True ) )
tableNameList = Group( delimitedList( tableName ) )

whereExpression = Forward()
and_ = Keyword("and", caseless=True)
or_ = Keyword("or", caseless=True)
in_ = Keyword("in", caseless=True)
# ARQ 2012-Feb-10: add "like" as an operator
like_ = Keyword("like", caseless=True)

E = CaselessLiteral("E")
# ARQ 2012-Feb-10: add "like" as a binop
binop = oneOf("= != < > >= <= eq ne lt le gt ge like", caseless=True)
arithSign = Word("+-",exact=1)
realNum = Combine( Optional(arithSign) + ( Word( nums ) + "." + Optional( Word(nums) ) |
                                                         ( "." + Word(nums) ) ) +
            Optional( E + Optional(arithSign) + Word(nums) ) )
intNum = Combine( Optional(arithSign) + Word( nums ) +
            Optional( E + Optional("+") + Word(nums) ) )

columnRval = realNum | intNum | quotedString | columnName # need to add support for alg expressions
whereCondition = Group(
    ( columnName + binop + columnRval ) |
    ( columnName + in_ + "(" + delimitedList( columnRval ) + ")" ) |
    ( columnName + in_ + "(" + selectStmt + ")" ) |
    ( "(" + whereExpression + ")" )
    )
whereExpression << whereCondition + ZeroOrMore( ( and_ | or_ ) + whereExpression )

# define the grammar
selectStmt << ( selectToken +
                   ( '*' | columnNameList ).setResultsName( "columns" ) +
                   fromToken +
                   tableNameList.setResultsName( "tables" ) +
                   Optional( Group( CaselessLiteral("where") + whereExpression ), "" ).setResultsName("where") )

simpleSQL = selectStmt

# ARQ 2012-Feb-10: define SQL Lite comment format, and ignore them
sqlLiteComment = "#" + restOfLine
simpleSQL.ignore( sqlLiteComment )

def parse_sql(str):
    try:
        return simpleSQL.parseString( str )
    except ParseException, err:
        print " "*err.loc + "^\n" + err.msg
        print err


def test( str ):
    print str,"->"
    try:
        tokens = simpleSQL.parseString( str )
        print "tokens = ", tokens
        print "tokens.columns =", tokens.columns
        print "tokens.tables =", tokens.tables
        print "tokens.where =", tokens.where
    except ParseException, err:
        print " "*err.loc + "^\n" + err.msg
        print err
    print

# test( "SELECT * from XYZZY, ABC" )
# test( "select * from SYS.XYZZY" )
# test( "Select A from Sys.dual" )
# test( "Select A,B,C from Sys.dual" )
# test( "Select A, B, C from Sys.dual" )
# test( "Select A, B, C from Sys.dual, Table2 " )
# test( "Xelect A, B, C from Sys.dual" )
# test( "Select A, B, C frox Sys.dual" )
# test( "Select" )
# test( "Select &&& frox Sys.dual" )
# test( "Select A from Sys.dual where a in ('RED','GREEN','BLUE')" )
# test( "Select A from Sys.dual where a in ('RED','GREEN','BLUE') and b in (10,20,30)" )
# test( "Select A,b from table1,table2 where table1.id eq table2.id and chrom = 1 and (a>2 or b <1)" )
# 
# """
# Test output:
# >pythonw -u simpleSQL.py
# SELECT * from XYZZY, ABC ->
# tokens = ['select', '*', 'from', ['XYZZY', 'ABC']]
# tokens.columns = *
# tokens.tables = ['XYZZY', 'ABC']
# 
# select * from SYS.XYZZY ->
# tokens = ['select', '*', 'from', ['SYS.XYZZY']]
# tokens.columns = *
# tokens.tables = ['SYS.XYZZY']
# 
# Select A from Sys.dual ->
# tokens = ['select', ['A'], 'from', ['SYS.DUAL']]
# tokens.columns = ['A']
# tokens.tables = ['SYS.DUAL']
# 
# Select A,B,C from Sys.dual ->
# tokens = ['select', ['A', 'B', 'C'], 'from', ['SYS.DUAL']]
# tokens.columns = ['A', 'B', 'C']
# tokens.tables = ['SYS.DUAL']
# 
# Select A, B, C from Sys.dual ->
# tokens = ['select', ['A', 'B', 'C'], 'from', ['SYS.DUAL']]
# tokens.columns = ['A', 'B', 'C']
# tokens.tables = ['SYS.DUAL']
# 
# Select A, B, C from Sys.dual, Table2 ->
# tokens = ['select', ['A', 'B', 'C'], 'from', ['SYS.DUAL', 'TABLE2']]
# tokens.columns = ['A', 'B', 'C']
# tokens.tables = ['SYS.DUAL', 'TABLE2']
# 
# Xelect A, B, C from Sys.dual ->
# ^
# Expected 'select'
# Expected 'select' (0), (1,1)
# 
# Select A, B, C frox Sys.dual ->
# ^
# Expected 'from'
# Expected 'from' (15), (1,16)
# 
# Select ->
# ^
# Expected '*'
# Expected '*' (6), (1,7)
# 
# Select &&& frox Sys.dual ->
# ^
# Expected '*'
# Expected '*' (7), (1,8)
# 
# >Exit code: 0
# """