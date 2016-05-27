import TableObjects

class CalibratorSource:

    def __init__(self, caltype, field_name, field_id, scan_ids):
        """
        
        """
        
        # Calibrator type
        self.caltype = caltype
        
        # Associated scan ids, field_id and field_name
        self.field_id = field_id
        self.scan_ids = scan_ids # list of strings
        self.field_name = field_name
        
        # Path extension (a string used for easy directory creation/access)
        # The files associated with this calibrator are saved in subdirectories
        # typically named as /flux_cal_name
        self.extend_dir = '/'+caltype+'_cal_'+self.field_name
        
        # String with the scan ids (casa format)
        self.scans = ','.join(self.scan_ids)
        
        # List of calibration tables derived from this calibrator
        self.derived_cal_tables = []
        
        
class TargetSource:

    def __init__(self, field_name, field_id, scan_ids):
        """
        
        """
        
        # Associated scan ids, field_id and field_name
        self.field_id = field_id
        self.scan_ids = scan_ids # list of strings
        self.field_name = field_name
        
        # Path extension (a string used for easy directory creation/access)
        # The files associated with this calibrator are saved in subdirectories
        # typically named as /flux_cal_name
        self.extend_dir = '/self_cal_'+self.field_name
        
        # String with the scan ids (casa format)
        self.scans = ','.join(self.scan_ids)
        
        # List of selc-cal tables derived from this target
        self.self_cal_gp_tables = []
        self.self_cal_ga_tables = []
