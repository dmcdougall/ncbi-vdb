/*===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 */

#include <klib/rc.h>
#include <klib/out.h>
#include <klib/text.h>
#include <klib/log.h>
#include <kproc/thread.h>

#include <xfs/xfs.h>
#include <xfs/tree.h>
#include "xfs-priv.h"
#include "schwarzschraube.h"
#include "xlog.h"

#include <stdlib.h>     /* using malloc() */
#include <string.h>     /* using memset() */
#include <unistd.h>     /* using STDOUT_FILENO, STDERR_FILENO */

/*  Some platform dependent headers
 */

#include "operations.h"

#ifndef STDIN_FILENO
#define STDIN_FILENO    0   /* Standard input.  */
#define STDOUT_FILENO   1   /* Standard output.  */
#define STDERR_FILENO   2   /* Standard error output.  */
#endif /* STDIN_FILENO */

/*
 *  Virtuhai table and it's methods
 */
static rc_t XFS_FUSE_init_v1 ( struct XFSControl * self );
static rc_t XFS_FUSE_destroy_v1 ( struct XFSControl * self );
static rc_t XFS_FUSE_mount_v1 ( struct XFSControl * self );
static rc_t XFS_FUSE_loop_v1 ( struct XFSControl * self );
static rc_t XFS_FUSE_unmount_v1 ( struct XFSControl * self);

static struct XFSControl_vt_v1 XFSControl_VT_V1 = {
    1,
    1,
    XFS_FUSE_init_v1,
    XFS_FUSE_destroy_v1,
    XFS_FUSE_mount_v1,
    XFS_FUSE_loop_v1,
    XFS_FUSE_unmount_v1
};

static struct fuse_operations TheFuseOperations;

/*  Control init.
 */
LIB_EXPORT
rc_t CC
XFSControlPlatformInit ( struct XFSControl * self )
{
    XFS_CAN ( self )

    self -> vt = ( union XFSControl_vt * ) & XFSControl_VT_V1;

    return 0;
}   /* XFSControlPlatformInit () */

/*  Overloadable versions
 */
rc_t
XFS_FUSE_init_v1 ( struct XFSControl * self )
{
    rc_t RCt = 0;

    XFS_CAN ( self )

    XFSLogMsg ( "XFS_FUSE_init()\n" );

    if ( self -> Control != NULL ) {

        return XFS_RC ( rcUnexpected );
    }

    if ( self -> Arguments == NULL ) {
        XFSLogErr ( XFS_RC ( rcUnexpected ), "XFS_FUSE_init(): arguments are empty\n" );
        return XFS_RC ( rcUnexpected );
    }

    if ( XFSControlGetLabel ( self ) == NULL ) {
        RCt = XFSControlSetLabel ( self, "FUSE" );
    }

    return RCt;
}   /* XFS_FUSE_init_v1 () */

rc_t
XFS_FUSE_destroy_v1 ( struct XFSControl * self )
{
    XFSLogMsg ( "XFS_FUSE_destroy()\n" );

    XFS_CAN ( self )

    return 0;
}   /* XFS_FUSE_destroy_v1 () */

rc_t
XFS_FUSE_mount_v1 ( struct XFSControl * self )
{
    rc_t RCt;

    struct fuse_args FuseArgs;
    struct fuse_chan * FuseChannel;
    struct fuse * FuseStruct;
    char * MountPoint;
    const char * LogPath;
    int Foreground;
    int Multithreaded;
    int Result;

    RCt = 0;
    FuseChannel = NULL;
    FuseStruct = NULL;
    MountPoint = NULL;
    LogPath = NULL;
    Foreground = true;
    Multithreaded = true;
    Result = 0;

    XFSLogMsg ( "XFS_FUSE_mount()\n" );

    XFS_CAN ( self )

    RCt = XFS_Private_InitOperations ( & TheFuseOperations );
    if ( RCt != 0 ) {
        XFSLogErr (
                XFS_RC ( rcFailed ),
                "Can not initialize operations\n"
                );
        return XFS_RC ( rcFailed ); 
    }

    memset ( & FuseArgs, 0, sizeof FuseArgs );
    Result = fuse_opt_add_arg ( & FuseArgs, XFSControlGetLabel ( self ) );
    Result = fuse_opt_add_arg ( & FuseArgs, XFSControlGetMountPoint ( self ) );
    LogPath = XFSControlGetLogFile ( self );

        /* Foreground */
    if ( XFSControlGetArg ( self, "-f" ) != NULL ) {
        Result = fuse_opt_add_arg ( & FuseArgs, "-f" );
    }

#if ! MAC
/* Options MAC does not know about */
    Result = fuse_opt_add_arg ( & FuseArgs, "-o" );
    Result = fuse_opt_add_arg( & FuseArgs, "big_writes" );
#endif /* MAC */

    if ( Result != 0 ) {
        XFSLogErr ( XFS_RC ( rcFailed ), "Can not mount\n" );
        return XFS_RC ( rcFailed ); 
    }

        /*  Adding the mountpoint
         */
    Result = fuse_parse_cmdline (
                    & FuseArgs,
                    & MountPoint,
                    & Multithreaded,
                    & Foreground
                );
    if ( Result != 0 ) {
        XFSLogErr ( XFS_RC ( rcFailed ), "Can not parse arguments\n" );
        return XFS_RC ( rcFailed ); 
    }

XFSLogDbg ( "Mnt = %s\nMlt = %d\nFrg = %d\n", MountPoint, Multithreaded, Foreground );


    FuseChannel = fuse_mount ( MountPoint, & FuseArgs );
    if ( FuseChannel == NULL ) {
        fuse_opt_free_args ( & FuseArgs );

        XFSLogErr ( XFS_RC ( rcFailed ), "Can not mount\n" );
        return XFS_RC ( rcFailed ); 
    }
    self -> ControlAux = FuseChannel;

        /*  Note passing TreeDepot as private data to fuse_context
         */
    FuseStruct = fuse_new (
                        FuseChannel,
                        & FuseArgs,
                        & TheFuseOperations,
                        sizeof ( struct fuse_operations ),
                        self -> TreeDepot
                    );
    if ( FuseStruct == NULL ) {
        fuse_unmount ( MountPoint, FuseChannel );

        fuse_opt_free_args ( & FuseArgs );

        XFSLogErr ( XFS_RC ( rcFailed ), "Can not fuse_new\n" );
        return XFS_RC ( rcFailed ); 
    }
    self -> Control = FuseStruct;

    Result = fuse_daemonize ( Foreground );
    if ( Result == -1 ) {
        fuse_unmount ( MountPoint, FuseChannel );

        fuse_opt_free_args ( & FuseArgs );

        XFSLogErr ( XFS_RC ( rcFailed ), "Can not daemonize\n" );
        return XFS_RC ( rcFailed ); 
    }

        /*  Here we are setting the log file
         */
    if ( LogPath != NULL ) {
        RCt = XFSLogInit ( LogPath );
        if ( RCt != 0 ) {
            fuse_unmount ( MountPoint, FuseChannel );

            fuse_opt_free_args ( & FuseArgs );

            XFSLogErr ( XFS_RC ( rcFailed ), "Can not set log file\n" );
            return RCt;
        }
    }

    Result = fuse_set_signal_handlers ( fuse_get_session ( FuseStruct ) );
    if ( Result == -1 ) {
        fuse_unmount ( MountPoint, FuseChannel );

        fuse_opt_free_args ( & FuseArgs );

        XFSLogErr (
                XFS_RC ( rcFailed ),
                "Can not set signal handlers\n"
                );
        return XFS_RC ( rcFailed ); 
    }

/*
    if ( Multithreaded != 0 ) {
        Result = fuse_loop_mt ( FuseStruct );
    }
    else {
        Result = fuse_loop ( FuseStruct );
    }

    fuse_remove_signal_handlers ( fuse_get_session ( FuseStruct ) );
*/

    fuse_opt_free_args ( & FuseArgs );

    XFSLogDbg ( " [XFS_FUSE_mount_v1] [%d] [%d]\n", __LINE__, Result );

    return Result == 0 ? 0 : XFS_RC ( rcInvalid );
}   /* XFS_FUSE_mount_v1 () */

rc_t
XFS_FUSE_loop_v1( struct XFSControl * self )
{
    rc_t RCt;
    struct fuse * FuseStruct;
    int Result;

    RCt = 0;
    FuseStruct = NULL;
    Result = 0;

    XFSLogMsg ( "XFS_FUSE_loop()\n" );

    if ( self == NULL ) {
        XFSLogErr ( 
                XFS_RC ( rcInvalid ),
                "XFS_FUSE_loop(): empty control passed\n"
                );
        RCt = 1;
    }
    else {
        FuseStruct = (struct fuse * ) self -> Control;

        Result = fuse_loop_mt ( FuseStruct );

        RCt = Result == 0 ? 0 : 1;
    }

    return RCt;
}   /* XFS_FUSE_loop_v1 () */

rc_t
XFS_FUSE_unmount_v1 ( struct XFSControl * self )
{
    struct fuse * FuseStruct = NULL;

    FuseStruct = (struct fuse * ) self -> Control;

    XFSLogDbg ( "XFS_FUSE_unmount()\n" );

    if ( self -> Control != NULL ) {

#if ! MAC
        fuse_exit ( FuseStruct );
#endif /* MAC */

XFSLogDbg ( "|o|fuse_unmount()\n" );
        fuse_unmount (
                    XFSControlGetMountPoint ( self ),
                    self -> ControlAux
                    );

XFSLogDbg ( "|o|fuse_remove_signal_handlers()\n" );
        fuse_remove_signal_handlers ( fuse_get_session ( FuseStruct ) );

XFSLogDbg ( "|o|fuse_destroy()\n" );
        fuse_destroy ( FuseStruct );

        self -> Control = NULL;
        self -> ControlAux = NULL;

XFSLogDbg ( "|o|exiting fuse()\n" );

        XFSLogDestroy ( );
    }
    else {
        XFSLogErr ( XFS_RC ( rcNull ), "XFS_FUSE_unmount(): empty control passed\n" );
    }

    return 0;
}   /* XFS_FUSE_unmount_v1 () */
