#!/usr/bin/env ruby
require 'optparse'
require 'statistics2'
require 'fileutils'

class Klass
  WORKDIR = File.expand_path("..", __dir__)
  OUTD = File.join(WORKDIR, "data")
  LOGD = File.join(WORKDIR, "logs")

  # Configurable, check the header of the previous step to adjust these values
  SAMPLE_IDX  = 11 # index of the first sample's GT column (0-based)
  SAMPLE_SB = 3 # index of the SB field within the sample-specific columns (0-based, relative to SAMPLE.GT)
  SAMPLE_COLS = 5 # number of columns with sample-sepcfic fields (GT, DP, AD, SB, GQ)
  attr_reader :opts

  def initialize(opts)
    @opts = opts
  end

  ADDITIONAL_HEADER = %w(Sample GT DP AD SB GQ)
  PHRED_MAX = Float::MIN_10_EXP * -10

  def phred_score(p)
    case
    when p > 1
      0
    when p > 0
      (Math.log10(p) * -10).round
    when p <= 0
      PHRED_MAX
    end
  end

  def chisq_test_phred(a, b, c, d)
    total = a + b + c + d
    row1_total = a + b
    row2_total = c + d
    col1_total = a + c
    col2_total = b + d
    e_a = (row1_total * col1_total).to_f / total
    e_b = (row1_total * col2_total).to_f / total
    e_c = (row2_total * col1_total).to_f / total
    e_d = (row2_total * col2_total).to_f / total
    chisq =
      ((a - e_a)**2 / e_a) + ((b - e_b)**2 / e_b) + ((c - e_c)**2 / e_c) + ((d - e_d)**2 / e_d)
    phred_score(Statistics2.chi2_x(1, chisq))
  end

  def boolianize(head:, name:, cols:)
    idx = head.index(name)
    if cols[idx] == "."
        cols[idx] = "false"
    else
        cols[idx] = "true"
    end
    cols
  end
  
  def simplify(fin:, fout:)
    header = nil
    snames = nil
    fin.each_line.with_index do |row, nrow|
      if nrow.zero?
        header = row.chomp.sub(/^\#/, '').split("\t")
        fout.puts header[0 .. (SAMPLE_IDX - 1)].concat(ADDITIONAL_HEADER).join("\t")
        snames = header[SAMPLE_IDX .. -1].each_slice(SAMPLE_COLS).map{|x|x[0][4 .. -5]}
        next
      end
      cols = row.chomp.split("\t")
      cols = boolianize(head: header, cols: cols, name:"RepeatMasker")
      cols = boolianize(head: header, cols: cols, name:"WindowMaskerSdust")
      cols = boolianize(head: header, cols: cols, name:"SegDup")
      cols = boolianize(head: header, cols: cols, name:"CpGisland")
      cols = boolianize(head: header, cols: cols, name:"SimpleRepeat")
      snames.each_with_index do |sname, idx|
        cols2 = Array.new
        cols2 << sname.dup
        gt_idx = SAMPLE_IDX + (SAMPLE_COLS * idx)
        # next if (cols[gt_idx] == "0/0" || cols[gt_idx].include?("."))
        (0 .. (SAMPLE_COLS - 1)).each do |i|
          if i == SAMPLE_SB
            cols2 << chisq_test_phred(*cols[gt_idx + i].split(",").map{|x|Float(x)})
            ### .to_s + "|" + cols[gt_idx + i].split(",").map{|x|x.to_i}.to_s
          else
            cols2 << cols[gt_idx + i]
          end
        end
        fout.puts (cols[0 .. (SAMPLE_IDX - 1)] + cols2).join("\t")
      end
    end
  end
  
  def run
    FileUtils.mkdir_p(LOGD)
    Dir["#{OUTD}/*neighbor2.annot.txt"].sort.each do |a2txt|
      out_base = File.basename(a2txt).sub(/\.annot\.txt$/, '.annot.sw.txt')
      out = File.join(OUTD, out_base)
      temp_out = File.join(OUTD, "__TEMP__#{out_base}")
      log = File.join(LOGD, "#{out_base}.log")
      open(log, 'w') do |flog|
        flog.puts("input=#{a2txt}")
        flog.puts("output=#{out}")
      end
      open(a2txt, 'r') do |fin|
        open(temp_out, 'w') do |fout|
          simplify(fin: fin, fout: fout)
        end
      end
      FileUtils.mv(temp_out, out)
    end
  end
end

if $0 == __FILE__
  RubyVM::YJIT.enable if defined? RubyVM::YJIT.enable
  opts = ARGV.getopts('h')
  Klass.new(opts.map{|k,v|[k.to_sym, v]}.to_h).run
end
