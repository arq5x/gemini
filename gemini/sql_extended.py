# select_parser.py
# Copyright 2010, Paul McGuire
#
# a simple SELECT statement parser, taken from SQLite's SELECT statement
# definition at http://www.sqlite.org/lang_select.html
#
from pyparsing import *

LPAR,RPAR,COMMA = map(Suppress,"(),")
select_stmt = Forward().setName("select statement")

# keywords
(UNION, ALL, AND, INTERSECT, EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, 
 CROSS, LEFT, OUTER, JOIN, AS, INDEXED, NOT, SELECT, DISTINCT, FROM, WHERE, GROUP, BY,
 HAVING, ORDER, BY, LIMIT, OFFSET) =  map(CaselessKeyword, """UNION, ALL, AND, INTERSECT, 
 EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, CROSS, LEFT, OUTER, JOIN, AS, INDEXED, NOT, SELECT, 
 DISTINCT, FROM, WHERE, GROUP, BY, HAVING, ORDER, BY, LIMIT, OFFSET""".replace(",","").split())
(CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, END, CASE, WHEN, THEN, EXISTS,
 COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, CURRENT_TIME, CURRENT_DATE, 
 CURRENT_TIMESTAMP) = map(CaselessKeyword, """CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, 
 END, CASE, WHEN, THEN, EXISTS, COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, 
 CURRENT_TIME, CURRENT_DATE, CURRENT_TIMESTAMP""".replace(",","").split())
keyword = MatchFirst((UNION, ALL, INTERSECT, EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, 
 CROSS, LEFT, OUTER, JOIN, AS, INDEXED, NOT, SELECT, DISTINCT, FROM, WHERE, GROUP, BY,
 HAVING, ORDER, BY, LIMIT, OFFSET, CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, END, CASE, WHEN, THEN, EXISTS,
 COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, CURRENT_TIME, CURRENT_DATE, 
 CURRENT_TIMESTAMP))
 
identifier = ~keyword + Word(alphas, alphanums+"_")
collation_name = identifier.copy()
column_name = identifier.copy()
column_alias = identifier.copy()
table_name = identifier.copy()
table_alias = identifier.copy()
index_name = identifier.copy()
function_name = identifier.copy()
parameter_name = identifier.copy()
database_name = identifier.copy()

# expression
expr = Forward().setName("expression")

integer = Regex(r"[+-]?\d+")
numeric_literal = Regex(r"\d+(\.\d*)?([eE][+-]?\d+)?")
string_literal = QuotedString("'")
blob_literal = Combine(oneOf("x X") + "'" + Word(hexnums) + "'")
literal_value = ( numeric_literal | string_literal | blob_literal |
    NULL | CURRENT_TIME | CURRENT_DATE | CURRENT_TIMESTAMP )
bind_parameter = (
    Word("?",nums) |
    Combine(oneOf(": @ $") + parameter_name)
    )
type_name = oneOf("TEXT REAL INTEGER BLOB NULL")

expr_term = (
    CAST + LPAR + expr + AS + type_name + RPAR |
    EXISTS + LPAR + select_stmt + RPAR |
    function_name + LPAR + Optional(delimitedList(expr)) + RPAR |
    literal_value |
    bind_parameter |
    identifier
    )

UNARY,BINARY,TERNARY=1,2,3
expr << operatorPrecedence(expr_term,
    [
    (oneOf('- + ~') | NOT, UNARY, opAssoc.LEFT),
    ('||', BINARY, opAssoc.LEFT),
    (oneOf('* / %'), BINARY, opAssoc.LEFT),
    (oneOf('+ -'), BINARY, opAssoc.LEFT),
    (oneOf('<< >> & |'), BINARY, opAssoc.LEFT),
    (oneOf('< <= > >='), BINARY, opAssoc.LEFT),
    (oneOf('= == != <>') | IS | IN | LIKE | GLOB | MATCH | REGEXP, BINARY, opAssoc.LEFT),
    ('||', BINARY, opAssoc.LEFT),
    ((BETWEEN,AND), TERNARY, opAssoc.LEFT),
    ])

compound_operator = (UNION + Optional(ALL) | INTERSECT | EXCEPT)

ordering_term = expr + Optional(COLLATE + collation_name) + Optional(ASC | DESC)

join_constraint = Optional(ON + expr | USING + LPAR + Group(delimitedList(column_name)) + RPAR)

join_op = COMMA | (Optional(NATURAL) + Optional(INNER | CROSS | LEFT + OUTER | LEFT | OUTER) + JOIN)

join_source = Forward()
single_source = ( (Group(database_name("database") + "." + table_name("table")) | table_name("table")) + 
                    Optional(Optional(AS) + table_alias("table_alias")) +
                    Optional(INDEXED + BY + index_name("name") | NOT + INDEXED)("index") | 
                  (LPAR + select_stmt + RPAR + Optional(Optional(AS) + table_alias)) | 
                  (LPAR + join_source + RPAR) )

join_source << single_source + ZeroOrMore(join_op + single_source + join_constraint)

result_column = "*" | table_name + "." + "*" | (expr + Optional(Optional(AS) + column_alias))

# ARQ 2012-Feb-10: add "like" as a binop
whereExpression = Forward()
and_ = Keyword("and", caseless=True)
or_ = Keyword("or", caseless=True)
in_ = Keyword("in", caseless=True)
like_ = Keyword("like", caseless=True)
E = CaselessLiteral("E")
binop = oneOf("= != < > >= <= eq ne lt le gt ge like", caseless=True)
arithSign = Word("+-",exact=1)
realNum = Combine( Optional(arithSign) + ( Word( nums ) + "." + Optional( Word(nums) ) |
                                                         ( "." + Word(nums) ) ) +
            Optional( E + Optional(arithSign) + Word(nums) ) )
intNum = Combine( Optional(arithSign) + Word( nums ) +
            Optional( E + Optional("+") + Word(nums) ) )

columnRval = realNum | intNum | quotedString | result_column # need to add support for alg expressions
ident = Word( alphas, alphanums + "_$" + ".$" ).setName("identifier")
columnName = delimitedList( ident, ".", combine=True ).setParseAction( upcaseTokens )
fromToken = Keyword("from", caseless=True)
tableName = delimitedList( ident, ".", combine=True ).setParseAction( upcaseTokens )
tableNameList = Group( delimitedList( tableName ) )
whereCondition = Group(
    ( columnName + binop + columnRval ) |
    ( columnName + in_ + "(" + delimitedList( columnRval ) + ")" ) |
    ( columnName + in_ + "(" + select_stmt + ")" ) |
    ( "(" + whereExpression + ")" )
    )
whereExpression << whereCondition + ZeroOrMore( ( and_ | or_ ) + whereExpression )

# Combine the grammar rules to define how a SQL statement should be parsed.
select_core = (SELECT + Optional(DISTINCT | ALL) + 
                ( '*' | Group(delimitedList(columnName)))("select") +
                fromToken +
                tableNameList.setResultsName("tables") +
                Optional(WHERE + Group(delimitedList(whereExpression))("where")))

# handle any order by, group by, having, or limit clauses
select_stmt << (select_core + ZeroOrMore(compound_operator + select_core) +
                Optional(ORDER + BY + Group(delimitedList(ordering_term))("orderby")) +
                Optional(GROUP + BY + Group(delimitedList(ordering_term))("groupby")) +
                Optional(HAVING + expr("having")) +
                Optional(LIMIT + (integer + OFFSET + integer | integer + COMMA + integer)))

simpleSQL = select_stmt

# ARQ 2012-Feb-10: define SQL Lite comment format, and ignore them
sqlLiteComment = "#" + restOfLine
simpleSQL.ignore( sqlLiteComment )

def parse_sql(str):
    try:
        return simpleSQL.parseString( str )
    except ParseException, err:
        print " "*err.loc + "^\n" + err.msg
        print err


# tests = """\
#     select * from variants
#     select * from variants where ref = alt
#     select * from xyzzy where z > 100 order by zz, yy
#     select * from xyzzy where z > 100 and x < 100 order by zz
#     select distinct chrom, start from variants where x > 100 and y < 200 order by zz
#     select distinct chrom, start from variants where x > 100 and y < 200 order by yy group by zz,xx having ww > 10
#     select *, start from variants where x > 100 and y < 200 order by zz
#     select * from xyzzy""".splitlines()
# for t in tests:
#     print t
#     try:
#         print select_stmt.parseString(t).dump()
#     except ParseException, pe:
#         print pe.msg
#     print
