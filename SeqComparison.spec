/*
A KBase module: SeqComparison
This sample module contains one small method - filter_contigs.
*/

module SeqComparison {
    /*
        Run DNAdiff.

        workspace_name - the name of the workspace for input/output
        input_genomeset_ref - optional input reference to genome set
        input_genome_refs - optional input list of references to genome objects
        output_report_name - the name of the output report

        n_break - break matches at n_break or more consecutive non-ACGTs

        @optional input_genomeset
        @optional input_genome_names
    */
    typedef structure {
        string workspace_name;
        string input_genomeset;
        list<string> input_genome_names;
        string output_report_name;
    } DNAdiffParams;

    typedef structure {
        string report_name;
        string report_ref;
    } DNAdiffOutput;

    funcdef run_mugsy(DNAdiffParams params) returns (DNAdiffOutput output)
        authentication required;

};
