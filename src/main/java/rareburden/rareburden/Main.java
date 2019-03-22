/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rareburden.rareburden;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.docopt.Docopt;

/**
 * Compute the gene burden score for each individual, variant list for genes were 
 * @author wallace
 */
public class Main {
    private static BufferedReader group_file_reader = null;
    private static VCFFileReader vcfreader = null;
    private static List<String> samplelist = null;
    private static Set<String>  sampleset = null;
    private static final StringJoiner sjoiner = new StringJoiner(":");
    private static final DecimalFormat FORMAT = new DecimalFormat("#.####E0");
    private static double max_maf = 1;
    
    private static final String DOC =
                "Compute gene burden score for each individual.\n"
                + "\n"
                + "Usage:\n"
                + "  RareBurdenScore -g file -v bgzfile --test txt [-s file] [--max-maf double] [-t cpus]\n"
                + "  RareBurdenScore (-h | --help)\n"
                + "  RareBurdenScore --version\n"
                + "\n"
                + "---------------------------\n"
                + "1. Output results to stdout.\n"
                + "2. Imupte the vcf missing value as ref allele.\n"
                + "---------------------------\n"
                + "\n"
                + "Options:\n"
                + "  -g file       Gene group file, with option inlcude vairant weight!\n"
                + "                mark-id as chr:pos:ref:alt, 4-elements no weight, 5 with weight.\n"
                + "                eg.: group_name 20:1110696:A:G,T:2 20:1234590:G:GTC:2\n"
                + "                     group_name 20:14370:G:A 20:1234590:G:GTC\n"
                + "  -v bgzfile    Bgziped and tabix indexed vcf file.\n"
                + "  -s file       Sample list file, one line per sample. !\n"
                + "                Load all samples from vcf if this option is off.\n"
                + "  --test txt    Specify the ways to compute burden score for each individual.\n"
                + "                b.collapse: count the number of rare(alt) alleles for each individual.\n"
                + "                b.collapse.weight: sum(alt_count*variant_weight).\n"
                + "  --max-maf double Only keep variant with maf <= max-maf.\n"
                + "  -t cpus       Number of cpus for computing.\n"
                + "  -h --help     Show this screen.\n"
                + "  --version     Show version.\n"
                + "\n";

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Map<String, Object> opts =
        new Docopt(DOC).withVersion("0.1").parse(args);

        if(opts.get("-t") != null){   
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
        }
        
        try {
//            String gfile = "/Users/wallace/NetBeansProjects/RareBurden/Test/test.group.txt";
            String gfile = (String) opts.get("-g");
            group_file_reader = new BufferedReader(new FileReader(gfile));
            
//            String[] namearray = {"NA00003","NA00002"};
            // load sample list from file.
            if(opts.get("-s") != null){   
                String[] namearray = new BufferedReader(new FileReader((String) opts.get("-s")))
                        .lines()
                        .filter(s->s.trim().length()>0)
                        .toArray(String[]::new);

                samplelist = Arrays.asList(namearray);
                sampleset = new HashSet<>(samplelist);
            }
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }

        
//        String vcfString = "/Users/wallace/NetBeansProjects/RareBurden/Test/test.indels.missing.vcf.gz";
        String vcfString = (String) opts.get("-v");
        File VCF = new File(vcfString);
        File VCF_TBI = new File(vcfString + ".tbi");
        
//        System.out.println("RareBurdenScore.Main.main()");
        
        vcfreader  = new VCFFileReader(VCF, VCF_TBI, true);
        VCFHeader vh = vcfreader.getFileHeader(); 
        // if no sample subset, use all samples from vcf file.
        if (samplelist == null) {
            samplelist = vh.getGenotypeSamples();
            sampleset = new HashSet<>(samplelist);
        }
        
        Set vcf_sampleSet = new HashSet(vh.getGenotypeSamples());
        long error_sample_count = samplelist.stream()
                .filter(s->!vcf_sampleSet.contains(s))
                .mapToInt(s->{
                    System.err.println("ERROR: can not file sample in VCF file, sample name: " + s);
                    return 1;
                })
                .count();
        if (error_sample_count>0) {System.exit(-1);}

        
        if(opts.get("--max-maf") != null){ max_maf = Double.parseDouble((String)opts.get("--max-maf"));}
                
        
        //output header.
        StringJoiner sj = new StringJoiner("\t");
        sj.add("#CHROM").add("BEGIN").add("END").add("MARKER_ID")
          .add("NUM_ALL_VARS").add("NUM_PASS_VARS").add("NUM_SING_VARS");
        samplelist.stream().forEach(s->sj.add(s));
        System.out.println(sj.toString());
        
        String TEST = (String)opts.get("--test");
        if (TEST.equalsIgnoreCase("b.collapse")) {
            run(false);
        }else if (TEST.equalsIgnoreCase("b.collapse.weight")) {
            run(true);
        }else{
            System.err.println("ERROR: Please set a proper value for --test");
            System.exit(-1);
        }
//        run();
    }
    
    /**
     * Compute the collapse score the each individual in each group.
     * @param with_weight : true: collapse.weight, false: collapse without weight (eaqual weight of 1 for each snp.).
     */
    private static void run(boolean with_weight) {
        //iterate gene by gene
        group_file_reader.lines()
           .filter(l -> l.trim().length()>0)
           .parallel()
           .forEach(group->{
               collapseScore(decodeGroupString(group),with_weight); //count the number of rare alleles.
           });
    }
    
    /**
     * Compute the aggregate rare_burden for each gene. Output results to stdout.
     * @param decodedGroup
     * @param with_weight : weight rare variants by the weight in group file. The fifth element.
     */
    private static void collapseScore(String[][] decodedGroup, boolean with_weight){
        String[] gname = decodedGroup[0];
//        System.err.println("Group ID: " + gname[0]);
//        AtomicDouble[] ID_SUM_SCORES = new AtomicDouble[samplelist.size()];
//        for (int i = 0; i < ID_SUM_SCORES.length; i++) {
//            ID_SUM_SCORES[i] = new AtomicDouble(0);
//        }

        //**** don't use parallel version, very strange behavior oberved.
        double[] ID_SUM_SCORES = new double[samplelist.size()]; Arrays.fill(ID_SUM_SCORES, 0.0);
        final int [] NUM_SINGLETON = new int[1]; NUM_SINGLETON[0] = 0;  //recode the number of singleton for this group.
        final int [] NUM_GOOD = new int[1]; NUM_GOOD[0] = 0; //number of sites in the scoring (passed QC).
        // NUM_ALL_VARS: is the number of variant in grouping file, but also found in VCF.
        final int [] NUM_ALL_VARS = new int[1]; NUM_ALL_VARS[0] = 0;
        
        final List<Integer> varpos = new ArrayList<>(decodedGroup.length-1);
        final Set<String> varidSet = new HashSet<>(decodedGroup.length -1);
        final Map<String, Double> snp_weight_map = new HashMap<>(decodedGroup.length-1);
        
        Arrays.stream(decodedGroup)
            .skip(1)
            .forEach(variant->{
                // varid as chr:pos:ref:alt 
                // replace "," by ":", to deal with multi-alelic allels issues.
                String varid = (variant[0] + ":" + variant[1] + ":" + variant[2] + ":" + variant[3]).replace(",",":").toUpperCase(Locale.ENGLISH);
                int pos = Integer.parseInt(variant[1]);
                varpos.add(pos);
                varidSet.add(varid);
                
                if (with_weight) {
                    snp_weight_map.put(varid, Double.parseDouble(variant[4]));
                }
            }); 
        
        String chr = decodedGroup[1][0];
        int min_pos = Collections.min(varpos);
        int max_pos = Collections.max(varpos);
        if (varidSet.size()>0) {
            //load a region for faster IO.
            //if the variant is indel with possible to load multiple sites. if the query pos is overlapped with.
            //make it as paralel version here.
            vcfreader.query(chr, min_pos, max_pos).stream()
                    //** don't inivoke parallel here, strange behavior observed.
//                    .parallel()
                    .forEach(s->{

//                System.err.println("------------------------------------------------------------------------");
                //loading the variants, and compute the aggregation of scores.
                //subset to the samples we want to load.
                VariantContext vc = s.subContextFromSamples(sampleset);
                StringJoiner sj = new StringJoiner(":");
                
//                for (Object string : vc.getAlternateAlleles().stream().map(a->a.getBaseString()).toArray()) {
                //*** We got all the alt alleles from all samples, in order to dealwith subsample all missing problem.
                //*** If all individual with missing allele here, alt will also be empty by htsjdk.
                for (Object string : s.getAlternateAlleles().stream().map(a->a.getBaseString()).toArray()) {
                    sj.add((String) string);
                }
                // deal with the problem, missing alt allele. alt allele coded as '.'
                // *** If all individual with missing allele here, alt will also be empty by htsjdk.
//                System.err.println("alt-alleles: " + sj.toString());
                if (sj.length()==0) {sj.add(".");}


                String loadVarID = vc.getContig() + ":" + vc.getStart() + ":" + vc.getReference().getBaseString() 
                        + ":" + sj.toString();
                loadVarID = loadVarID.toUpperCase(Locale.ENGLISH);
//                System.err.println("loadVarID: " + loadVarID);


                // match id.     
                if (!varidSet.contains(loadVarID)) { return; }
                NUM_ALL_VARS[0]+=1;
                
                //compute and filter by max-maf: 
//                System.err.println("TOTAL_CALLED\t" + vc.getCalledChrCount()); //number number of called alleles.
//                System.err.println("TOTAL_CALLED_REF\t"+vc.getCalledChrCount(vc.getReference())); //number of called ref. alleles.
    //                    int total_count = vc.getHomRefCount() + vc.getHetCount() + vc.getHomVarCount();
    //                    if (total_count == 0) {return;}
    //                    double af = (vc.getHetCount()/2.0 + vc.getHomVarCount())/total_count;
                int num_called_alleles = vc.getCalledChrCount();
                int num_called_ref_allels = vc.getCalledChrCount(vc.getReference());
                if (num_called_alleles-num_called_ref_allels == 1) {
                    NUM_SINGLETON[0] +=1;
//                      NUM_SINGLETON[0].addAndGet(1);
                }

                double af =  num_called_ref_allels * 1.0 / num_called_alleles;
                double maf = Math.min(af, 1-af);
//                System.err.println("MAF: " + maf);
                //by defual max_maf if 1, therefore there are no maf filterring.
                if (maf==0 || maf > max_maf ) {return;}

                //Compute the scores.
                double snp_weight = 1.0;
                if (with_weight) {
                    snp_weight = snp_weight_map.get(loadVarID);
                }
                
                //cumulation of scores for each individual.
                Iterable<Genotype> genotypes = vc.getGenotypesOrderedBy(samplelist);
                int i = 0; 
                for (Genotype g : genotypes) {
                    double score = (2-g.countAllele(vc.getReference())) * snp_weight;
                    if (g.isNoCall()) { score = 0;} // missing set the weight as ref_alleles.
                    ID_SUM_SCORES[i++] += score;
//                    System.err.println("score: " + score);
//                      ID_SUM_SCORES[i++].addAndGet(score);
                }

                NUM_GOOD[0] += 1;
//                NUM_GOOD[0].addAndGet(1);
            });
        }
        
//        sj.add("#CHROM").add("BEGIN").add("END").add("MARKER_ID")
//          .add("NUM_ALL_VARS").add("NUM_PASS_VARS").add("NUM_SING_VARS");
        
        StringJoiner sj = new StringJoiner("\t");
        sj.add(chr).add(Integer.toString(min_pos)).add(Integer.toString(max_pos));
        sj.add(gname[0]);
        // add the number of vaiants in group file, passed QC, singleton.
        // NUM_ALL_VARS: is the number of variant in grouping file, but also found in VCF.
//        sj.add(Integer.toString(decodedGroup.length-1));
        sj.add(Integer.toString(NUM_ALL_VARS[0]));
        sj.add(Integer.toString(NUM_GOOD[0]));
        sj.add(Integer.toString(NUM_SINGLETON[0]));
        // output the sum score for each indiviauls.
        Arrays.stream(ID_SUM_SCORES)
                .forEach(s->{sj.add(FORMAT.format(s));});
        System.out.println(sj.toString());
    }
    
    /**
     * decode group file. 
     * @param gString
     * @return 2D array, first line as group name. 2-n line, element attributes of a group.
     */
    private static String[][] decodeGroupString(String gString){
        String[] garr = gString.split("\\s+");
        String[][] results = new String[garr.length][];
        String[] gname = new String[1]; gname[0] = garr[0];
        results[0] = gname;
        for (int i = 1; i < garr.length; i++) {
            results[i] = garr[i].split(":");
        }
//        for (String[] result : results) {
//            System.out.println(Arrays.toString(result));
//        }
        return results;
    }
    
}
