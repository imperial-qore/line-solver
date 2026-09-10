classdef LineStatus
    % LINESTATUS One line of output that is REWRITTEN rather than repeated.
    %
    % An iterative solver reports its progress once per iteration. Printed as
    % one line each, a run of a few hundred iterations scrolls everything else
    % out of the terminal to say the same six numbers over and over. LINESTATUS
    % keeps that report on a SINGLE row and updates it in place:
    %
    %   LineStatus.set('Iter %2d. Analyze time: %.3fs.', it, t);
    %   LineStatus.append(' MaxIterErr=%.3e', err);   % same row, more fields
    %   LineStatus.close();                           % end the row, once
    %
    % The row is rewound with BACKSPACES, not with a carriage return: '\b' is
    % the idiom that works both in the MATLAB desktop Command Window and under
    % `matlab -batch`, whereas '\r' is handled inconsistently between them. A
    % shorter row is padded with blanks so the tail of a longer predecessor
    % cannot survive underneath it.
    %
    % AN INTERLEAVED PRINT ENDS THE ROW, it does not corrupt it. Backspaces
    % walk back over whatever is actually on the terminal and cannot cross a
    % newline, so a warning raised mid-iteration used to land inside the row
    % and leave the NEXT rewind short, splicing two iterations into one line:
    %
    %   Iter 51. ... Runtime: 6.177s. MaxIterErr=1.08e-14 (tol=5.0000Iter 52. ...
    %
    % LINE_PRINTF therefore calls CLOSE before it writes anything, which is why
    % this class emits through its own raw path rather than through line_printf
    % -- otherwise every update would close the row it was drawing. A warning
    % now gets its own line and the next iteration opens a fresh row. Inside
    % SolverLN the interleaving is rare to begin with, because the layer
    % solvers are silenced (see SolverLN.setSolver).
    %
    % CLOSE is idempotent and emits the newline the row never had.
    %
    % A row containing a newline cannot be rewound, so SET and APPEND strip
    % them. At VerboseLevel.SILENT every entry point is a no-op and no state is
    % kept, so the bookkeeping cannot drift against what was actually written.
    %
    % See also: line_printf, LineConsole, VerboseLevel
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)

        function st = state(newst)
            % ST = STATE(NEWST) reads, and optionally replaces, the row state
            %
            % TWO FIELDS, AND THEY ARE NOT THE SAME STRING. `text` is the row's
            % LOGICAL content, which APPEND extends; `width` is how many
            % characters are actually on the terminal, which is what the next
            % rewind has to walk back over. Keeping only the on-screen string
            % and appending to that made every append start after the blanks
            % that padded the previous row, so the row grew by its own padding
            % on every iteration and ran to hundreds of thousands of columns.
            persistent ROWSTATE
            if isempty(ROWSTATE)
                ROWSTATE = struct('text', '', 'width', 0);
            end
            if nargin > 0
                ROWSTATE = newst;
            end
            st = ROWSTATE;
        end

        function tf = isOpen()
            % TF = ISOPEN() true while a status row is being rewritten
            tf = LineStatus.state().width > 0;
        end

        function reset()
            % RESET() forgets the open row WITHOUT emitting anything
            %
            % For an interrupted solve, where the row is already lost: closing
            % it would backspace over whatever the error printed instead.
            LineStatus.state(struct('text', '', 'width', 0));
        end

        function set(fmt, varargin)
            % SET(FMT, ...) replaces the row with this text
            LineStatus.render(LineStatus.clean(sprintf(fmt, varargin{:})));
        end

        function append(fmt, varargin)
            % APPEND(FMT, ...) extends the row and rewrites it
            %
            % The row is assembled by more than one caller -- EnsembleSolver
            % ITERATE lays down the timings and SolverLN CONVERGED adds the
            % iteration error -- so appending has to re-render the whole row
            % rather than print a fragment after it.
            LineStatus.render([LineStatus.state().text, ...
                LineStatus.clean(sprintf(fmt, varargin{:}))]);
        end

        function close()
            % CLOSE() ends the row with a newline and forgets it
            %
            % The blanks that padded a SHRINKING row stay on the finished line.
            % They cannot be taken back: a backspace moves the cursor, it does
            % not erase, so rewinding over the pad and then emitting the
            % newline changes only where the cursor was, not what the line
            % holds. They are invisible on a terminal and only surface if the
            % line is copied. The way to avoid them is to keep one-off notices
            % OFF the row -- see how SolverLN.converged prints its averaging
            % notice on a line of its own -- so that the row never shrinks.
            if LineStatus.isOpen()
                LineStatus.reset();
                LineStatus.emit(sprintf('\n'));
            end
        end

    end

    methods (Static, Access = private)

        function emit(s)
            % EMIT(S) writes S to the output stream WITHOUT going through
            % line_printf
            %
            % line_printf now closes an open row before it prints anything
            % (that is how a warning gets its own line instead of landing in
            % the middle of this one), so a row that wrote through it would
            % close itself on every update. This is the raw path; SILENT is
            % checked by the callers, which is also where the bookkeeping
            % lives.
            fprintf(GlobalConstants.StdOut, '%s', s);
        end

        function s = clean(s)
            % S = CLEAN(S) a row is one line: backspaces cannot cross a newline
            s = strrep(strrep(s, sprintf('\n'), ' '), sprintf('\r'), ' ');
        end

        function render(txt)
            % RENDER(TXT) rewinds the previous row and writes this one
            if GlobalConstants.Verbose == VerboseLevel.SILENT
                return % nothing was written, so nothing is tracked
            end
            prev = LineStatus.state();
            % PAD TO THE PREVIOUS WIDTH so a row that got shorter does not
            % leave the tail of the longer one standing after it. The padding
            % is written but NOT remembered as text: only its width is, so an
            % APPEND extends the row rather than the blanks after it.
            shown = txt;
            if numel(shown) < prev.width
                shown = [shown, blanks(prev.width - numel(shown))];
            end
            LineStatus.emit([repmat(sprintf('\b'), 1, prev.width), shown]);
            LineStatus.state(struct('text', txt, 'width', numel(shown)));
        end

    end
end
