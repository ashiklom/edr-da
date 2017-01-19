load_template <- function(example_xml_path = '~/Projects/pecan/outputs/PEcAn_77000000004/pecan.xml',
                          template_path = 'pecan_template.xml') {
    file.copy(example_xml_path, template_path, overwrite = TRUE)
}

new_settings <- function(outdir, workflowid, css, pss, site, pft_list,
                      template_path = 'pecan_template.xml') {
    stopifnot(is.character(outdir),
              is.character(workflowid),
              is.character(css),
              is.character(pss),
              is.character(site),
              is.character(pft_list),
              file.exists(template_path))

    stopifnot(length(pft_list > 0),
              all(names(pft_list) == 'pft'))

    template_settings <- XML::xmlToList(XML::xmlParse(file = template_path))

    newsettings <- template_settings

    # Set output directory
    newsettings$outdir <- outdir

    # Set workflowID
    newsettings$workflow$id <- workflowid

    # Set inputs
    newsettings$inputs$css$id <- css
    newsettings$inputs$pss$id <- pss
    newsettings$inputs$site$id <- site

    # Set PFTs
    newsettings$pfts <- pft_list
    # May also need to set PFT number...?

    return(new_settings)
}

write_new_settings <- function(new_settings, outfile) {
    xml <- PEcAn.utils::listToXml(x = new_settings, tag = 'pecan')
    XML::saveXML(doc = xml, file = outfile)
    return(xml)
}
