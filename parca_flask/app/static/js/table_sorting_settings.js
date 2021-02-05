
$(function() {

    $(".tablesorter").tablesorter({
      
  
      widthFixed : true,
      widgets: ["zebra", "filter"],
      ignoreCase: false,
      widgetOptions : {
  
        filter_childRows : false,
        filter_childByColumn : false,
        filter_childWithSibs : true,
        filter_columnFilters : true,
        filter_columnAnyMatch: true,
        filter_cellFilter : '',
        filter_cssFilter : '', // or []
        filter_defaultFilter : {},
        filter_excludeFilter : {},
        filter_external : '',
        filter_filteredRow : 'filtered',
        filter_formatter : null,
        filter_functions : null,
        filter_hideEmpty : false,
        filter_hideFilters : false,
        filter_ignoreCase : true,
        filter_liveSearch : true,
        filter_matchType : { 'input': 'exact', 'select': 'exact' },
        filter_onlyAvail : 'filter-onlyAvail',
        filter_placeholder : { search : '', select : '' },
        filter_reset : 'button.reset',
        filter_resetOnEsc : true,
        filter_saveFilters : true,
        filter_searchDelay : 300,
        filter_searchFiltered: true,
        filter_selectSource  : null,
        filter_serversideFiltering : false,
        filter_startsWith : false,
        filter_useParsedData : false,
        filter_defaultAttrib : 'data-value',
        filter_selectSourceSeparator : '|'
  
      }
  
    });
  
    // Clear stored filters - added v2.25.6
    $('.resetsaved').click(function(){
      $('.tablesorter').trigger('filterResetSaved');
  
      // show quick popup to indicate something happened
      var $message = $('<span class="results"> Reset</span>').insertAfter(this);
      setTimeout(function(){
        $message.remove();
      }, 500);
      return false;
    });
    $('button[data-filter-column]').click(function(){
    
      var filters = [],
        $t = $(this),
        col = $t.data('filter-column'), // zero-based index
        txt = $t.data('filter-text') || $t.text(); // text to add to filter
  
      filters[col] = txt;
    
      $.tablesorter.setFilters( $('.tablesorter'), filters, true ); // new v2.9
  
     
      var filters = $('table.tablesorter').find('input.tablesorter-filter');
      filters.val(''); // clear all filters
      filters.eq(col).val(txt).trigger('search', false);
  
      var columns = [];
      columns[5] = '2?%'; // or define the array this way [ '', '', '', '', '', '2?%' ]
      $('table').trigger('search', [ columns ]);
  
      return false;
    });
  
  });
   $('.header').on("click",function(){
     $(this).nextUntil('tr.header').toggle();
     $(this).closest('table').find('tbody').toggle();
     $(this).closest('table').find('tfoot').toggle();
     $(this).find('span').text(function(_, value){return value=='-'?'+':'-'});
  });
  
  
  