#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "export.h"
#include "idl_aacgm.h"

/*
 * Define message codes and their corresponding printf(3) format
 * strings. Note that message codes start at zero and each one is
 * one less that the previous one. Codes must be monotonic and
 * contiguous. Must match the corresponding entries in IDLtoC.h
 */
static IDL_MSG_DEF msg_arr[] =
{
  {  "idl_aacgm_ERROR",		"%NError: %s." },
  {  "idl_aacgm_NOSTRINGARRAY",    "%NString arrays not allowed %s"},
 };

/*
 * The load function fills in this message block handle with the
 * opaque handle to the message block used for this module. The other
 * routines can then use it to throw errors from this block.
 */
IDL_MSG_BLOCK msg_block;

int IDL_Load(void)
{
/*define the message block*/
  if (!(msg_block = IDL_MessageDefineBlock("AACGM", ARRLEN(msg_arr),
	   msg_arr))) {
	return IDL_FALSE;
  }
/*Call the startup function to add the routines to IDL.*/
  if (!idl_aacgm_startup()){
    IDL_MessageFromBlock(msg_block, idl_aacgm_ERROR, 
		IDL_MSG_RET,"Unable to initialize AACGM");
  }
  
  return IDL_TRUE;
}
