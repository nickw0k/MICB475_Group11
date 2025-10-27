# MICB475: Data Science Research in Microbiology 2025W Group 11

## Table of Contents 
  * [Meeting 1](#meeting-1)
  * [Meeting 2 Agenda](#meeting-2-Agenda)
  * [Meeting 2](#meeting-2)
  * [Meeting 3 Agenda](#meeting-3-Agenda)
  * [Meeting 3](#meeting-3)
  * [Meeting 4 Agenda](#meeting-4-Agenda)
  * [Meeting 4](#meeting-4)
  * [Meeting 5 Agenda](#meeting-5-Agenda)


    
## Meeting 1 
Date: Monday, Sept. 29 (2:00-2:45)

Minutes prepared by Nicholas Kwok.

People present: Evelyn, Hans, Jason, Luna, Thomas, Nicholas 

Notetaker: Nicholas Kwok

Expectations: 
- Write agenda every meeting. Can ask Hans over the weekend if we have any questions
- We need a notetaker every meeting (can rotate).
- Add Hans into the Github next week
- Has to be a novel research question. Have to finalize by end of this week if we are thinking of a new one. 

Notes: 
- We want a novel research question, something that hasn't been done before. 
- Use AI to find  datasets specific conditions?
- Parkinsons is saturated
- Multiple Sclerosis is also saturated?
- Alcohol Consumption dataset? Has not been widely used before, possible choice 
- Smoking and vaping dataset
- Write down possible questions, find niche questions online
- Pick something you're interested in!
- Has to be centered around the microbiome
- Bacteria mediated COVID immunity? (Nick)
- Gut microbiome Sepsis (Jason) 
- Smoking (Thomas) 
- Infant (Luna)
- Thursday, groupchat discussion after finding 10-15 datasets on our topic, then deciding on which topic we want
- ask in group if we have questions

Admin:
- Moved meetings to tuesdays for Meeting 2 and 3 (7th Oct, 14th Oct 2:00pm-2:45pm)
- Finish data 2 weeks before presentation (first meeting talk about story, talk about slides, second meeting make edits)


## Meeting 2 Agenda

Prepared by Nicholas Kwok

Date: Tuesday Oct 7 
Time: 2:00 – 2:45 PM 
Location: LSI 1510 
Notetaker: Luna

1. Administrative Updates
  - Confirm all members have access to the GitHub repo and uploaded “Possible Datasets” PDF. (Jason?)
  - Add Hans to the repo (verify permissions).
  - Review timeline: Finalize dataset choice this week → draft research question by Oct 14 meeting.

2. Dataset Evaluation 
  - Objective: Narrow from topic options (COVID, Infant Feeding, Smoking/Vaping, Sepsis) to one theme based on data quality and novelty.
  
    Discussion Prompts:
      - Evelyn’s feedback identified three COVID studies with usable metadata (Wei et al. 2023; Serrato et al. 2023; Newsome et al. 2021).
      - Compare COVID vs Infant Feeding vs Sepsis vs Smoking datasets for:
        Metadata completeness (accession IDs, n size, geographic range)
        Microbiome site (oral, gut, saliva, fecal)
        Novelty / previous analyses (saturation of literature)
        Feasibility of downloading and processing (16S vs metagenomic)
            
  - *Final Goal: Vote for top two datasets to pilot in QIIME2.*
  
3. Research Question Brainstorm 

Each member pitches 1–2 potential questions based on the shortlisted dataset(s).
Check that questions fit the “novel and analyzable” criteria.

4. Division of Work
  - Based on chosen theme (COVID/smoking/infant/Sepsis), Assign tasks for remaining dataset groups:

    Proposed division of work:
    - (If there are more than 1 explorable themes, 2 members analyze metadata + come up w/ novel research questions of each theme.)
    - (If there is only 1 explorable theme, 2 members try to come up with novel research questions, 2 analyze the metadata.
  
    Tasks before Meeting 3:
    - Confirm data availability & metadata quality.
    - Import one dataset into QIIME2 and generate a basic summary (ASV table, metadata linking).
    - Start a shared Google Doc for research-question drafting.
    
5. Next Steps & Deadlines 

  - Finalize chosen dataset by Oct 10 (Fri) for approval by Evelyn.
  - Write and upload a concise research question + hypothesis to GitHub by Oct 14 (Tues).
  - Next meeting (Oct 14) → review initial analyses in QIIME2. 



## Meeting 2
Date: Tuesday, Oct 7th (2:00-2:45)

Notetaker: Luna Kurokawa

People present: Evelyn, Hans, Jason, Luna, Thomas, Nicholas 

- Decide on which topic to do today
- Proposal is in two weeks so we want to ideally start thinking about the research question today

- 3 Covid supplemental papers was accepted - Evelyn does not agree with of one of the paper's analysis but the data itself was okay. Only diversity metrics were done
- We can do paired data using 2 of the covid datasets but it might be challenging due to the study location on the body of each paper (oral vs. gut). We don't have to combine it but anlayse it parallel (US and Mexico). (pre-analysis in qiime is separate) If one turns out to be a better dataset than the other we can just focus on one. Start with US dataset and Mexico dataset and if we have time contiue on to the China dataset. 
- Mexico looks saliva samples and has better annotiations 
- Jordan and Russia covid dataset do not have controls 
- Continuing on with the fecal analysis from the covid paper is anlso an option
- The big picture is to find micorbial signature of COVID patients in oral vs fecal microbiomes.

- By next week:
    - Commit on dataset 
    - Make a formal research question
    - What we want to explore in this topic
         - abundance? function of community? etc?
    
    - Start on qiime2 analysis
    - Hans will help us build a outline for the proposal
    - Next week we can start dicussing about what to analyse. 
    
## Meeting 3 Agenda
Created by Thomas Hesselbo
October 14th, 2025:  2:00pm - 2:45pm
LSI 1510
Notetaker: Thomas Hesselbo

1. Make the final decision on which dataset(s) to use
    - Decide if we are going to use the COVID, Infant nutrition, or Smoking dataset(s).
        - This was not spoken about during the last discussion as COVID dataset discussion dominated the conversation.
          If COVID, decide which of the 3 COVID data sets to use. I.E. Jordan + Mexico, if needed + China

2. Determine area of analysis within the datasets
    - Selection of parameters and approach to analysis
        - This will rely on analysis of previous papers to avoid covering the same topic (I.E., American and Mexico research question comprehensiveness).
3. Final Formulation of initial research question (It may change later on depending on the dataset results and so on).
    - Presentation of research questions made by group members
    - Decide if research question will be predictive or descriptive
        - (I.E., can x microbiome be used to determine likelihood of y disease vs x and y microbiomes have distinctive compositions based on z).
    - Largely dependent on the area of analysis and quality of data.
      
4. Proposal outline
   - Overall Timeline
   - Division of work
      - Citation
      - Denoising tool
      - Sample size + sequencing depth (pre and post denoising)
      - Max read length
      - Read length trends
      - Read truncation
      - Read quality plot
      - Sample selection
      - Rarefaction depth
      - Rarefaction curve
    
        
## Meeting 3

Participant Thomas, Luna, Jason, Nicholas, Hans

1. Dataset
   - Settle on COVID data sets
   - For the sub datasets, unsure:
2. Research question:
   - Usage of AI tool "undermine"
   - Predictive vs descriptive : a matter of preference, model training (predictive)
     - Using one dataset to predict data in another (testing vs training dataset) may cause problems if microbiomes are too different
   - Broad question: multiple aims for answering the question
3. Proposal timeline:
   - Need to decide on a research question ASAP as by next monday is too short of a timeframe before     proposal is due.

 4. Pre- pre meeting 4
   - Create code for rarefaction analysis before pre pre meeting
   - Two potential research questions per person on COVID dataset

 5. Pre meeting 4
    - Final proposed research question
    - Introduction (potentially, to be discussed during pre-pre meeting)
    - Sub dataset decision


Notes on COVID microbiome paper from Nicholas:
   In Newsome et al. and most early COVID microbiome papers, researchers only described taxonomic shifts (changes in abundance of genera/families).
They did not extend these findings to functional potential (what those microbes might be doing).
If you use the same 16S rRNA data but add a new functional annotation layer, for example:
Use PICRUSt2 to infer KEGG pathways from your USA cohort (Newsome dataset), or
Map predicted enzyme functions (MetaCyc pathways) and test how they differ between infected/recovered/control,
then you’re generating novel biological insight from existing data.
That’s legitimate novelty because:
The functional inference wasn’t done in the original paper.
You’re testing different biological questions (e.g., “Does recovery correlate with restoration of butyrate biosynthesis pathways?”).

## Meeting 4 Agenda

Prepared by Jason Jonathan

Date: Monday, Oct 20
Time: 2:00 – 2:45 PM 
Location: LSI 1416

**1. Finalize main research question for the project**
   - Discussed and decided on one prospective question to bring to Hans: Do age and sex-specific differences in the commensal microbiota influence the abundance and functional roles of microbial marker species that differentiate severe from mild COVID-19 patients?
     - One issue with this question: QIIME2 analysis beta-diversity metric showed that there no significant differences between male and female microbiome (p-value > 0.05)
     - Have back-up research questions ready in case the question above is not good enough:
       - Does functional redundancy among indicator taxa differ between infection and recovery states?
       - Are specific immune-related metabolic pathways (tryptophan, butyrate, bile acid metabolism) reliable functional indicators of recovery?

**2. Determine division of work for the proposal**
   - Decide on the parts of the proposal for each member to work on or if multiple members will on the same parts together
   - Decide project timeline for the proposal
   - Create a workflow overview chart for the proposal
  
**3. Decide on a time to meet and discuss outside of regular team meetings before proposal deadline**

**4. Next steps and deadlines**
   - Start writing team proposal and aim to finish by Saturday, Oct 25 latest to allow time for the team to review proposal together
   - Team proposal deadline Sunday, Oct 26 11:59PM
    
**5. Questions to ask Hans:**
   - What p-value are we willing to accept for sex-specific differences?
   - How do we define the age groups (ages range from teens to 80+)?

## Meeting 4

Date: Monday, Oct 20 (2:00 - 2:45 PM) in LSI 1416

Notetaker: Jason Jonathan

Participants: Jason, Nicholas, Luna, Thomas, Hans

**1. Main research question:** Do age and sex-specific differences in the commensal microbiota influence the abundance and functional roles of microbial marker species that differentiate severe from mild COVID-19 patients?
   - Advice from Hans:
     - Too specific with "commensal"; microbiome sequencing data may show both commensal and pathogenic bacteria
     - Reframe "sex-specific differences in the commensal microbiota"; this can mean sex-specific differences in the microbiota itself rather than the sex of the patients
     - Remove "microbial marker species"; we don't know what are the microbial marker species yet, instead we will find out through our analysis
     - Hans suggested to NOT do age for our analysis, however if we still wanted to we can define specific age groups we can work with (these need to be justified properly)
       - If we decide to do age for our question, there will be more variables we need to analyze which would make it more difficult; Hans said with sex and mild/severe COVID variables there would be plenty enough for our analysis
       - Most of the patients in the dataset are middle-aged, so if we are doing age for our analysis it's likely to be 40+
   - ***FINAL RESEARCH QUESTION (REFRAMED BY HANS): Are there sex-specific differences in microbial taxonomic and functional diversity in patients with severe or mild COVID?***

**2. Project proposal (due Sunday, Oct 26 11:59PM):**
   - Start thinking about aims for the proposal before the weekend
     - Aim 1 could be alpha and beta diversity analysis
     - Aim 2 could be core microbiome analysis, indicator species analysis, or DESeq analysis
     - Aim 3 could be functional diversity analysis with Picrust2 or differentail abundance analysis
     - Choose what analysis to do for each aim based on how the analysis we do can help us answer the research question, make sure each aim helps us answer the research question in one way or another
     - For writing proposal, one column per aim and each row lists the specific step we would do to answer the corresponding aim (can refer to Proposal_Example1.pdf on MICB 475 Canvas page module)
   - Start finding appropriate reference papers for the proposal
     - Make sure that the reference papers chosen are relevant to the question and helps us set the background in the introduction as well as directs our analysis
     - Hans suggested a minimum of 20 reference papers
   - Use Proposal_Example1.pdf from Canvas page as a guideline for writing the proposal (can be found in Project 2 Assignment Documents module)

**3. Schedule a meeting with Hans before the submission deadline**
   - Nick said maybe Thursday or Friday?

**4. Schedule a discord call between us to discuss proposal writing?**
   - Nick and Luna will do Introduction & Background?
   - Jason and Thomas will do the code and analysis (Dataset Overview)?
   - Will come up and do aims together as well as Gantt chart and Project Timeline

## Meeting 5 Agenda 

Prepared by Nicholas Kwok

Date: Monday, Oct 27
Time: 2:00 – 2:45 PM 
Location: LSI 1416

**1. Quick recap**
- Everyone gives a quick update on progress 
- Hans can review current progress on proposal doc

**2. Confirm 3 aims** 

Aim 1 — Assess whether sex influences overall microbial community structure within mild and severe COVID-19 groups

- Discuss alpha & beta diversity metrics (Shannon, Faith’s PD, weighted/unweighted UniFrac PCoA).
- Confirm statistical comparisons (PERMANOVA, pairwise tests).
- Update on QIIME2 diversity core metrics and R results by Thomas/Jason

Aim 2 — Identify sex-specific microbial taxa associated with disease severity

- Outline workflow for indicator species analysis (ISA): input = phyloseq object from mpt_final.
- Confirm thresholds (min abundance/prevalence).
- Decide whether to add DESeq2 for validation.
- Assign coding responsibilities for ISA and visualization (Venn/heatmaps). (Luna?)

Will ask Hans for advice?

Aim 3: Evaluate whether predicted functional pathways differ by sex and disease severity

- Plan PICRUSt2 analysis (ASV table + taxonomy -> predict KEGG orthology/ pathways/ Enzymes)
- Nick will handle pathway differential analysis and figures (PCoA + bar plots)

**3. Discuss deadlines**

- Since deadline got extended:

- Proposed goals:
1. Meeting on Sunday evening (5pm ish) 
2. Finished own parts by 6pm Saturday Nov 1st








