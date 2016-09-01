#!/usr/bin/env ruby

# Adapted from Stack Overflow answer http://superuser.com/a/738951

# "This script takes a single multiple-page pdf file as first argument,
# and outputs a latex file that multiplexes those pages
# Requires that pdfinfo is installed"

require 'optparse'

# Note to self: Ruby's option parser is weird
options = {}
rest = OptionParser.new do |opt|
    opt.banner = \
"This script takes a multiple-page pdf file, and outputs the latex code that multiplexes those pages.
Output can be piped to pdflatex to produce the multiplexed pdf.
Requires that pdfinfo is installed.

USAGE: pdftile.rb [options] PDF file [| pdflatex]"

    opt.on('--per-page INT', Numeric, 'Number of images per page') do |i| 
        options[:per_page] = i
    end
    opt.on('--per-row INT', Numeric, 'Number of images in a row') do |i| 
        options[:per_row] = i
    end
    opt.on('PDF file to multiplex') do |f|
    end
end.parse!
options[:filename] = rest[0]

raise OptionParser::MissingArgument.new("You must set a number of plots per page with --per-page") if options[:per_page].nil?
raise OptionParser::MissingArgument.new("You must set a number of plots per row with --per-row") if options[:per_row].nil?

nimg = (options[:per_page]).to_i
ncol = (options[:per_row]).to_i
if (nimg % ncol).zero?
    nrow = nimg / ncol
else
    nrow = (nimg / ncol) + 1
end

imgheight = 0.975 / nrow.to_f  # Slightly less than 1/nrow to allow for margins
imgwidth = 0.975 / ncol.to_f

latexhead = <<'EOF'
\documentclass{article}
\usepackage[pdftex]{graphicx}
\usepackage[margin=0.1in]{geometry}
\usepackage{pdfpages}
\begin{document}
EOF
latextail = <<'EOF'
\end{document}
EOF

pages = %x[pdfinfo #{options[:filename]}].split(/\n/).select{|x| x=~ /Pages:/}[0].split(/\s+/)[1].to_i
puts latexhead
s = (1..pages).each_slice(options[:per_page]).to_a
s.each do |a|
  puts "\\begin{figure}[!ht]"
  a.each do |p|
    puts "\\includegraphics[page=#{p},scale=0.2,height=#{'%.03f' % imgheight}\\textheight,width=#{'%.03f' % imgwidth}\\textwidth]{#{options[:filename]}}"
  end
  puts "\\end{figure}"
end
puts latextail