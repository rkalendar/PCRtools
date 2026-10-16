/* #####################
# 
#  Notice (do not remove or edit):
#  This JavaScript is created by Kaare Byberg.
#
#  BEGIN VALIDATE FORM  
#     
##################### */

/*
  USAGE:

  E.g.:
    <form action="#" method="post" onsubmit="if( validate('first_name', 'First Name','r',0)==false) return false;">

  * If r: Required field
  * If n: Must be a number
  * If 0: Must be 1> character
  * If another number, e.g. 5, it has to match the length exactly
  * If e, it's an e-mail address

  It only took 2 hours to figure out :-)

  --EspenA
  
  Important note: If any of the names in the validation call are incorrect, the user will be forwarded to 
  the URL  found in the action="" instead of being told that something is wrong.
*/

function validate() 
  {
  var i; 
  var thearguments = validate.arguments; 

  var errors = '';
  var errormessage = '';
  
  var either = '';
  var fieldgroup = '';
  var sensor = 1;

    for (i = 0; i < (thearguments.length - 3); i += 4) 
    { 
    var minimumsize = thearguments[i+3];
    var tests = thearguments[i+2];  
    var fieldname = thearguments[i+1];
    var field = document.getElementById(thearguments[i]);
    var contents = field.value;
    var size = contents.length;

    if (tests.indexOf('g') != -1) either += '\n\t' + fieldname;

    if (contents != "") 
      {
      if (tests.indexOf('g') != -1) sensor *= 0;

      if (tests.indexOf('e') != -1) 
        { 
          var mailaddress = /^(\S*)@{1}(\S*)$/;
          if ((contents.match(mailaddress)) == null) errors += '\n• ' + fieldname + ' must be an e-mail address.';
        } 
      else if (tests.indexOf('n') != -1) 
        {           
          if (isNaN(contents)) errors += '\n• ' + fieldname + ' must be a number.';      
        }

      if (size < minimumsize) errors += '\n• ' + fieldname + ' value is too short.'; 
      } 
    else if (tests.indexOf('r') != -1) errors += '\n• ' + fieldname + ' is a required field.'; 
    else if (tests.indexOf('g') != -1) sensor *= 1;
    }  
  if (sensor == 1 && either != '') fieldgroup += 'At least one of these fields must be filled in:' + either;

  if (errors) errormessage += 'The following fields are incorrect:' + errors;

  if (errormessage && fieldgroup) fieldgroup += '\n\n'
 
  if (errormessage || fieldgroup) {
    alert(fieldgroup + errormessage);
    return false;
  }
  else
    return true;
    
  //document.validateok = (errormessage == '' && fieldgroup == '');
  }

/* #####################
#                      
#  END VALIDATE FORM  
#                      
##################### */

