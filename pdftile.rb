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
    opt.on('--page-width FLOAT', Numeric, 'Output page width in mm [Optional - scanned from input pdf]') do |i|
        options[:page_width] = i
    end
    opt.on('--page-height FLOAT', Numeric, 'Output page height in mm [Optional - scanned from input pdf]') do |i|
        options[:page_height] = i
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

# Scan pdf file with pdfinfo to get number of pages and dimensions
v = %x[pdfinfo #{options[:filename]}].split(/\n/).select{|x| x=~ /Pages:/ or x=~ /Page size:/}.map { |x| x.split(/\s+/) }
pages = v[0][1].to_i
if options[:page_width].nil?
    original_width = v[1][2].to_f / 2.845  # Convert to mm (~2.845 pts per mm)
    options[:page_width] = ncol * original_width
end
if options[:page_height].nil?
    original_height = v[1][4].to_f / 2.845
    options[:page_height] = nrow * original_height
end


latexhead = <<-EOF
\\documentclass{article}
\\usepackage[pdftex]{graphicx}
\\usepackage[margin=0.1in,paperheight=#{'%.02f' % options[:page_height]}mm,paperwidth=#{'%.02f' % options[:page_width]}mm]{geometry}
\\usepackage{pdfpages}
\\begin{document}
EOF
latextail = <<'EOF'
\end{document}
EOF

puts latexhead
s = (1..pages).each_slice(options[:per_page]).to_a
s.each do |a|
  puts "\\begin{figure}[!ht]"
  a.each do |p|
    puts "\\includegraphics[page=#{p},height=#{'%.03f' % imgheight}\\textheight,width=#{'%.03f' % imgwidth}\\textwidth]{#{options[:filename]}}"
    puts "\\centering"
  end
  puts "\\end{figure}"
end
puts latextail
