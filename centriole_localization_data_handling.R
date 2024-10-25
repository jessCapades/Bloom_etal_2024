## centriole localization data handling



## stack bar chart

create_stacked_bar_chart = function (summed_data, selected_stage, fill_condition) {
  
  visualization = ggplot(filter(summed_data, stage == selected_stage),
                         mapping = aes(fill = localized_properly,
                                       x = embryo_type,
                                      y = Count))
  
  visualization + geom_bar(position="fill", stat="identity") +
    xlab("Cross Type") + ylab("Fraction of Total Embryos") +
    theme(text = element_text(size = 15)) +
    theme(panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
    theme(legend.text=element_text(size = 11)) + scale_fill_brewer()
}

