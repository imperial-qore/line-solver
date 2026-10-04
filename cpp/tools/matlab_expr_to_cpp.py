#!/usr/bin/env python3
"""Translate a MATLAB scalar arithmetic expression into C++, mechanically.

The kpctoolbox MAP(2) fitters (map_mmpp2.m, map_block.m) carry Maple-generated
closed forms tens of kilobytes long. They must not be retyped: a single wrong
sign in 50 KB of algebra produces a plausible MAP that matches no moment, and no
review catches it. This parses the expression and re-emits it, so the C++ is a
mechanical image of the MATLAB.

Grammar: numbers, identifiers, unary +/-, binary + - * / ^ (^ right-associative
and binding tighter than unary minus, as in MATLAB), parentheses, and calls.
Element-wise operators .* ./ .^ are accepted and emitted as the scalar ones,
which is what they mean on scalars.

Emission:
  x^(1/2)  -> num_sqrt(x)          (the port's exact-arithmetic square root)
  x^k      -> pw(x, k)             for an integer literal k
  x^y      -> num_pow(x, y)        otherwise
  numbers  -> num_traits<T>::from_int / from_rational, never a double literal,
              so the expression stays exact under Rational arithmetic.

Usage:  matlab_expr_to_cpp.py <file.m> <lhs-name> [<lhs-name> ...]
prints, for each assignment `<lhs> = <expr>;` found in the file, the C++ form.
"""
import re
import sys


class Parser(object):
    def __init__(self, s):
        self.s = s
        self.i = 0
        self.n = len(s)

    def ws(self):
        while self.i < self.n and self.s[self.i] in ' \t\n\r.':
            # a lone '.' only ever appears as part of .* ./ .^ here
            if self.s[self.i] == '.' and self.i + 1 < self.n and self.s[self.i + 1] not in '*/^':
                break
            if self.s[self.i] == '.' and self.i + 1 < self.n and self.s[self.i + 1] in '*/^':
                break
            self.i += 1

    def peek(self):
        self.ws()
        return self.s[self.i] if self.i < self.n else ''

    def eat(self, c):
        self.ws()
        if self.s[self.i:self.i + len(c)] == c:
            self.i += len(c)
            return True
        return False

    def parse(self):
        e = self.expr()
        self.ws()
        if self.i != self.n:
            raise SyntaxError('trailing input at %d: %r' % (self.i, self.s[self.i:self.i + 40]))
        return e

    def expr(self):
        e = self.term()
        while True:
            self.ws()
            if self.eat('+'):
                e = '(%s + %s)' % (e, self.term())
            elif self.i < self.n and self.s[self.i] == '-':
                self.i += 1
                e = '(%s - %s)' % (e, self.term())
            else:
                return e

    def term(self):
        e = self.unary()
        while True:
            self.ws()
            if self.eat('.*') or self.eat('*'):
                e = '(%s * %s)' % (e, self.unary())
            elif self.eat('./') or self.eat('/'):
                e = '(%s / %s)' % (e, self.unary())
            else:
                return e

    def unary(self):
        self.ws()
        if self.eat('-'):
            return '(-%s)' % self.unary()
        if self.eat('+'):
            return self.unary()
        return self.power()

    def power(self):
        base = self.atom()
        self.ws()
        if self.eat('.^') or self.eat('^'):
            exp_src_start = self.i
            e = self.unary()  # right-associative, and -x binds inside the exponent
            return self.emit_pow(base, e, self.s[exp_src_start:self.i].strip())
        return base

    def emit_pow(self, base, expcpp, expsrc):
        src = expsrc.replace(' ', '')
        if src in ('(1/2)', '1/2', '0.5'):
            return 'num_sqrt(%s)' % base
        m = re.fullmatch(r'\(?(\d+)\)?', src)
        if m:
            return 'pw(%s, %s)' % (base, m.group(1))
        m = re.fullmatch(r'\(?-(\d+)\)?', src)
        if m:
            return '(num_traits<T>::from_int(1) / pw(%s, %s))' % (base, m.group(1))
        return 'num_pow(%s, %s)' % (base, expcpp)

    def atom(self):
        self.ws()
        if self.eat('('):
            e = self.expr()
            if not self.eat(')'):
                raise SyntaxError('missing ) at %d' % self.i)
            return e
        m = re.match(r'\d+\.\d+([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?|\d+[eE][+-]?\d+', self.s[self.i:])
        if m:
            self.i += m.end()
            return 'num_traits<T>::from_double(%s)' % m.group(0)
        m = re.match(r'\d+', self.s[self.i:])
        if m:
            self.i += m.end()
            return 'num_traits<T>::from_int(%s)' % m.group(0)
        m = re.match(r'[A-Za-z_]\w*', self.s[self.i:])
        if m:
            self.i += m.end()
            name = m.group(0)
            self.ws()
            if self.i < self.n and self.s[self.i] == '(':
                self.i += 1
                args = []
                if self.peek() != ')':
                    args.append(self.expr())
                    while self.eat(','):
                        args.append(self.expr())
                if not self.eat(')'):
                    raise SyntaxError('missing ) in call at %d' % self.i)
                fn = {'sqrt': 'num_sqrt', 'abs': 'num_abs'}.get(name, name)
                return '%s(%s)' % (fn, ', '.join(args))
            return name
        raise SyntaxError('unexpected %r at %d' % (self.s[self.i:self.i + 20], self.i))


def assignments(text, names):
    """Yield (name, expression-source) for `name = ...;` at statement start."""
    out = {}
    for name in names:
        m = re.search(r'^\s*' + re.escape(name) + r'\s*=\s*(.*?);\s*$', text,
                      re.MULTILINE | re.DOTALL)
        if not m:
            raise KeyError('no assignment to %s' % name)
        # stop at the first `;` that is not inside parentheses
        src, depth = [], 0
        for ch in m.group(1):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            elif ch == ';' and depth == 0:
                break
            src.append(ch)
        out[name] = ''.join(src)
    return out


def cse(exprs, minlen=40):
    """Bind repeated subexpressions to temporaries.

    Maple output repeats the same radical and the same denominator dozens of
    times; emitting it literally makes the C++ both unreadable and O(repeats)
    to evaluate. Identical balanced-parenthesis substrings longer than `minlen`
    that occur more than once across all the expressions become `cse<k>`.
    Purely textual, so it cannot change the value: the same string denotes the
    same value in the same scope.
    """
    def subexprs(s):
        stack, out = [], {}
        for i, ch in enumerate(s):
            if ch == '(':
                stack.append(i)
            elif ch == ')' and stack:
                a = stack.pop()
                sub = s[a:i + 1]
                if len(sub) >= minlen:
                    out[sub] = out.get(sub, 0) + 1
        return out

    counts = {}
    for e in exprs.values():
        for sub, n in subexprs(e).items():
            counts[sub] = counts.get(sub, 0) + n

    # Longest first, so an outer repeat is bound before the inner ones it holds.
    cands = sorted([s for s, n in counts.items() if n > 1], key=len, reverse=True)
    binds = []
    for sub in cands:
        total = sum(e.count(sub) for e in exprs.values()) + sum(b[1].count(sub) for b in binds)
        if total < 2:
            continue  # already absorbed into an earlier, longer binding
        name = 'cse%d' % len(binds)
        binds.append((name, sub))
        for k in list(exprs):
            exprs[k] = exprs[k].replace(sub, name)
        for i in range(len(binds) - 1):
            binds[i] = (binds[i][0], binds[i][1].replace(sub, name))
    # A binding may reference a later one; order them so each is defined first.
    ordered, emitted = [], set()

    def emit(name, body):
        if name in emitted:
            return
        for other, obody in binds:
            if other != name and re.search(r'\b%s\b' % other, body):
                emit(other, obody)
        emitted.add(name)
        ordered.append((name, body))

    for name, body in binds:
        emit(name, body)
    return ordered, exprs


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    text = open(sys.argv[1]).read()
    text = re.sub(r'\.\.\.\s*\n', '', text)      # MATLAB line continuations
    text = re.sub(r'^\s*%.*$', '', text, flags=re.MULTILINE)
    exprs = {}
    for name, src in assignments(text, sys.argv[2:]).items():
        exprs[name] = Parser(src).parse()
    binds, exprs = cse(exprs)
    for name, body in binds:
        print('    const T %s = %s;' % (name, body))
    if binds:
        print()
    for name in sys.argv[2:]:
        print('    const T %s = %s;' % (name, exprs[name]))
        print()


if __name__ == '__main__':
    main()
