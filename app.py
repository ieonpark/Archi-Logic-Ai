import streamlit as st
import math
import folium
from folium.plugins import Draw
from streamlit_folium import st_folium
from geopy.geocoders import Nominatim
from pyproj import Geod
from shapely.geometry import mapping, Polygon, MultiPolygon
import pandas as pd

st.set_page_config(layout="wide", page_title="Archi Logic AI - Feasibility Analysis", initial_sidebar_state="expanded")

# ============ METADATA ============
AUTHOR = "IEON PARK"
CONTACT_EMAIL = "IEONPARK@HOTMAIL.COM"
COPYRIGHT_YEAR = 2025

# ============ UTILITIES & DATA ============
geod = Geod(ellps="WGS84")
geocoder = Nominatim(user_agent="archi_logic_app")

fire_lane_requirements = {
    "Georgia": {"min": 20, "code": "IBC 503.2.1 + Georgia Amendments", "notes": "26 ft required if serving buildings >30 ft height"},
    "North Carolina": {"min": 20, "code": "NC Building Code 503.2.1", "notes": "26 ft for aerial apparatus access required"},
    "Texas": {"min": 26, "code": "TAS 503.2.1 (Texas Accessibility Standards)", "notes": "TxDOT may require 30 ft for certain classifications"}
}

lihtc_qap_defaults = {
    "Georgia": {
        "parking_reduction_allowed": True,
        "max_parking_reduction_pct": 20,
        "qap_strictness": "Medium",
    },
    "North Carolina": {
        "parking_reduction_allowed": False,
        "max_parking_reduction_pct": 0,
        "qap_strictness": "High",
    },
    "Texas": {
        "parking_reduction_allowed": True,
        "max_parking_reduction_pct": 25,
        "qap_strictness": "Low",
    }
}

def geocode_address(address):
    try:
        loc = geocoder.geocode(address, timeout=10)
        if not loc:
            return None
        return {"lat": loc.latitude, "lon": loc.longitude, "display_name": loc.address}
    except Exception:
        return None

def geodesic_area_sqft(geom):
    try:
        if isinstance(geom, (Polygon, MultiPolygon)):
            geom_mapping = mapping(geom)
        elif isinstance(geom, dict):
            geom_mapping = geom
        else:
            return 0.0
        gtype = geom_mapping.get("type")
        coords = geom_mapping.get("coordinates")
        if not coords:
            return 0.0
        total_m2 = 0.0
        if gtype == "Polygon":
            for ring in coords:
                if len(ring) < 3:
                    continue
                lons, lats = zip(*ring)
                area, _ = geod.polygon_area_perimeter(lons, lats)
                total_m2 += abs(area)
        elif gtype == "MultiPolygon":
            for poly in coords:
                for ring in poly:
                    if len(ring) < 3:
                        continue
                    lons, lats = zip(*ring)
                    area, _ = geod.polygon_area_perimeter(lons, lats)
                    total_m2 += abs(area)
        return total_m2 * 10.76391041671
    except Exception:
        return 0.0

def normalize_folium_return(fret):
    features = []
    if not fret:
        return features
    if isinstance(fret, dict):
        if "all_drawings" in fret and fret["all_drawings"]:
            cand = fret["all_drawings"]
            if isinstance(cand, dict) and "features" in cand:
                return cand["features"]
        for key in ("geojson", "last_draw", "last_drawn", "last_object", "last_geojson"):
            if key in fret and fret[key]:
                val = fret[key]
                if isinstance(val, dict) and val.get("type") == "FeatureCollection" and val.get("features"):
                    return val["features"]
                if isinstance(val, dict) and val.get("type") == "Feature":
                    return [val]
                if isinstance(val, dict) and val.get("geometry"):
                    return [val]
        for key in ("objects", "features", "all_features"):
            if key in fret and isinstance(fret[key], list):
                return fret[key]
        for v in fret.values():
            if isinstance(v, dict) and v.get("geometry"):
                return [v]
    return features

def try_area_with_orders(coords):
    area1 = area2 = 0.0
    try:
        lons1, lats1 = zip(*coords)
        a1, _ = geod.polygon_area_perimeter(lons1, lats1)
        area1 = abs(a1)
    except Exception:
        area1 = 0.0
    try:
        swapped = [(pt[1], pt[0]) for pt in coords]
        lons2, lats2 = zip(*swapped)
        a2, _ = geod.polygon_area_perimeter(lons2, lats2)
        area2 = abs(a2)
    except Exception:
        area2 = 0.0
    if area2 > area1:
        return swapped, area2 * 10.76391041671
    else:
        return coords, area1 * 10.76391041671

def line_to_polygon_geom_with_order_corr(geom):
    try:
        gtype = geom.get("type")
        if gtype == "LineString":
            coords = geom.get("coordinates", [])
            if len(coords) < 3:
                return None
            if coords[0] != coords[-1]:
                coords = coords + [coords[0]]
            best_coords, _ = try_area_with_orders(coords)
            return {"type": "Polygon", "coordinates": [best_coords]}
        if gtype == "MultiLineString":
            mcoords = geom.get("coordinates", [])
            if not mcoords or len(mcoords[0]) < 3:
                return None
            coords = mcoords[0]
            if coords[0] != coords[-1]:
                coords = coords + [coords[0]]
            best_coords, _ = try_area_with_orders(coords)
            return {"type": "Polygon", "coordinates": [best_coords]}
        return geom
    except Exception:
        return None

def meters_per_degree(lat):
    lat_rad = math.radians(lat)
    meters_per_deg_lat = 111132.92 - 559.82 * math.cos(2 * lat_rad) + 1.175 * math.cos(4 * lat_rad)
    meters_per_deg_lon = 111412.84 * math.cos(lat_rad) - 93.5 * math.cos(3 * lat_rad)
    return meters_per_deg_lat, meters_per_deg_lon

def scale_geometry_to_target(geom, target_area_sf):
    try:
        current_area = geodesic_area_sqft(geom)
        if current_area <= 0:
            return geom
        factor = math.sqrt(float(target_area_sf) / float(current_area))
        coords = geom.get("coordinates")
        if not coords:
            return geom
        if geom["type"] == "Polygon":
            rings = coords
            flat = rings[0]
            lons = [p[0] for p in flat]
            lats = [p[1] for p in flat]
            centroid_lon = sum(lons) / len(lons)
            centroid_lat = sum(lats) / len(lats)
            mlat, mlon = meters_per_degree(centroid_lat)
            new_rings = []
            for ring in rings:
                new_ring = []
                for (lon, lat) in ring:
                    dx_m = (lon - centroid_lon) * mlon
                    dy_m = (lat - centroid_lat) * mlat
                    dx_m *= factor
                    dy_m *= factor
                    new_lon = centroid_lon + (dx_m / mlon)
                    new_lat = centroid_lat + (dy_m / mlat)
                    new_ring.append([new_lon, new_lat])
                new_rings.append(new_ring)
            return {"type": "Polygon", "coordinates": new_rings}
        elif geom["type"] == "MultiPolygon":
            new_polys = []
            for poly in coords:
                flat = poly[0]
                lons = [p[0] for p in flat]
                lats = [p[1] for p in flat]
                centroid_lon = sum(lons) / len(lons)
                centroid_lat = sum(lats) / len(lats)
                mlat, mlon = meters_per_degree(centroid_lat)
                new_poly = []
                for ring in poly:
                    new_ring = []
                    for (lon, lat) in ring:
                        dx_m = (lon - centroid_lon) * mlon
                        dy_m = (lat - centroid_lat) * mlat
                        dx_m *= factor
                        dy_m *= factor
                        new_lon = centroid_lon + (dx_m / mlon)
                        new_lat = centroid_lat + (dy_m / mlat)
                        new_ring.append([new_lon, new_lat])
                    new_poly.append(new_ring)
                new_polys.append(new_poly)
            return {"type": "MultiPolygon", "coordinates": new_polys}
        else:
            return geom
    except Exception:
        return geom

# ============ SESSION STATE ============
if "workspace_features" not in st.session_state:
    st.session_state["workspace_features"] = []
if "role_map" not in st.session_state:
    st.session_state["role_map"] = {}
if "map_center" not in st.session_state:
    st.session_state["map_center"] = (38.9072, -77.0369)
if "map_zoom" not in st.session_state:
    st.session_state["map_zoom"] = 17

# ============ SIDEBAR INPUTS (COMMON) ============
st.sidebar.header("Project Parameters")
state = st.sidebar.selectbox("State", ["Georgia", "North Carolina", "Texas"])
project_type = st.sidebar.selectbox("Project Type", ["Multifamily General", "Senior Housing (55+)", "Mixed-Income"])
funding_type = st.sidebar.selectbox(
    "Funding Type",
    ["General", "LIHTC 9% (Competitive)", "LIHTC 4% (Tax-Exempt Bond)", "Other"]
)

st.sidebar.subheader("Unit Mix")
studio_units = st.sidebar.number_input("Studio Units", min_value=0, value=0, step=1)
br1_units = st.sidebar.number_input("1 Bedroom Units", min_value=0, value=20, step=1)
br2_units = st.sidebar.number_input("2 Bedroom Units", min_value=0, value=25, step=1)
br3_units = st.sidebar.number_input("3+ Bedroom Units", min_value=0, value=7, step=1)
total_units = studio_units + br1_units + br2_units + br3_units

fha_min_accessible_default = max(1, math.ceil(total_units * 0.05)) if total_units > 0 else 0
st.sidebar.subheader("Accessibility")
accessible_units = st.sidebar.number_input("Accessible Units (Total)", min_value=0, max_value=total_units if total_units>0 else 0, value=fha_min_accessible_default, step=1)
type_a_units = st.sidebar.number_input("Type A (Mobility) Units", min_value=0, max_value=accessible_units if accessible_units>0 else 0, value=accessible_units if accessible_units>0 else 0, step=1)
type_b_units = max(0, accessible_units - type_a_units)
fha_min_accessible = fha_min_accessible_default
fha_compliant = accessible_units >= fha_min_accessible

st.sidebar.subheader("Site Constraints")
site_area_sf = st.sidebar.number_input("Total Site Area (SF)", min_value=0, value=50000, step=1000)
building_footprint_target = st.sidebar.number_input("Building Footprint (SF) - target", min_value=0, value=15000, step=500)
available_parking_input = st.sidebar.number_input("Available Parking Spaces (count)", min_value=0, value=45, step=1)
site_width = st.sidebar.number_input("Fire Lane Width (ft)", min_value=0, value=20, step=1)

st.sidebar.markdown("---")
st.sidebar.subheader("QAP Scoring Factors (expanded)")

# Tier 1: Location & Access
st.sidebar.write("Location & Access")
transit_access = st.sidebar.checkbox("Transit within 0.5 mile", value=False, help="+10 points")
opportunity_zone = st.sidebar.checkbox("Opportunity Zone", value=False, help="+15 points")
high_opportunity_area = st.sidebar.checkbox("High opportunity area (job growth, low poverty)", value=False, help="+15 points")

walk_score = st.sidebar.slider("Walk Score", 0, 100, 50, help="0-5 points based on score")

nearby_services = st.sidebar.number_input("Essential services (1 mile)", 0, 10, 3, help="2 points each, max 10")

school_proximity = st.sidebar.checkbox("Public school within 1 mile", value=False, help="+5 points")
job_centers = st.sidebar.checkbox("Major employers/job centers within 3 miles", value=False, help="+10 points")
grocery_store = st.sidebar.checkbox("Grocery store within 1 mile", value=False, help="+5 points")
medical_facilities = st.sidebar.checkbox("Medical facilities within 1 mile", value=False, help="+5 points")
parks_recreation = st.sidebar.checkbox("Parks/recreation within 0.5 mile", value=False, help="+3 points")

# Tier 2: Site Readiness
st.sidebar.write("Site Readiness")
site_control_documented = st.sidebar.checkbox("Site control documented", value=True, help="+5 points")
environmental_clear = st.sidebar.checkbox("Phase I ESA clear", value=True, help="+8 points")
utilities_available = st.sidebar.checkbox("All utilities available at property line", value=True, help="+3 points")
zoning_compliant = st.sidebar.checkbox("Site properly zoned (no variance needed)", value=True, help="+5 points")

# Tier 3: Site Development Features
st.sidebar.write("Site Development Features")
green_space = st.sidebar.checkbox("Green space/outdoor amenities (10%+ of site)", value=False, help="+5 points")
playground = st.sidebar.checkbox("Playground equipment", value=False, help="+3 points")
community_building = st.sidebar.checkbox("Community center/meeting space", value=False, help="+10 points")
bike_storage = st.sidebar.checkbox("Covered bike storage", value=False, help="+2 points")

# Tier 4: Sustainability
st.sidebar.write("Sustainability & Efficiency")
energy_certification = st.sidebar.selectbox(
    "Energy Certification Plan",
    ["None", "Energy Star", "LEED Certified", "LEED Silver", "LEED Gold+"],
    help="0/5/10/12/15 points"
)
solar_ready = st.sidebar.checkbox("Solar ready or solar installed", value=False, help="+8 points")
water_conservation = st.sidebar.checkbox("Water conservation measures", value=False, help="+3 points")
ev_charging = st.sidebar.checkbox("EV charging station(s)", value=False, help="+3 points")

# Tier 5: Community Support
st.sidebar.write("Community Engagement")
community_support = st.sidebar.checkbox("Community support letters obtained", value=False, help="+5 points")
local_government = st.sidebar.checkbox("Local government support/resolution", value=False, help="+8 points")
service_coordination = st.sidebar.checkbox("Service provider coordination agreements", value=False, help="+5 points")

# Tier 6: Population Served
st.sidebar.write("Target Population (if applicable)")
special_needs_commitment = st.sidebar.checkbox("Commitment to serve special needs (10%+ units)", value=False, help="+10 points")

st.sidebar.markdown("---")
st.sidebar.markdown("Map / Drawing Settings")
map_height = st.sidebar.number_input("Map Height (px)", min_value=300, value=600, step=50)
prefer_drawn_parking = st.sidebar.checkbox("Prefer drawn parking (by area) over sidebar parking count", value=True)

STANDARD_SPACE_SF = 180
ACCESSIBLE_SPACE_SF = 216
VAN_SPACE_SF = 288

# ============ CALCULATIONS (COMMON) ============
def calculate_parking_requirement(studio, br1, br2, br3, project_type, state):
    if project_type == "Senior Housing (55+)":
        ratio_studio = 0.5; ratio_1br = 0.5; ratio_2br = 0.75; ratio_3br = 1.0
    else:
        if state == "Georgia":
            ratio_studio = 1.0; ratio_1br = 1.0; ratio_2br = 1.5; ratio_3br = 2.0
        elif state == "North Carolina":
            ratio_studio = 1.0; ratio_1br = 1.25; ratio_2br = 1.75; ratio_3br = 2.0
        else:
            ratio_studio = 1.0; ratio_1br = 1.5; ratio_2br = 2.0; ratio_3br = 2.5
    required_raw = (studio * ratio_studio + br1 * ratio_1br + br2 * ratio_2br + br3 * ratio_3br)
    required = int(math.ceil(required_raw))
    if required <= 0:
        ada_required = 0
    elif required <= 25:
        ada_required = 1
    elif required <= 50:
        ada_required = 2
    elif required <= 75:
        ada_required = 3
    elif required <= 100:
        ada_required = 4
    elif required <= 150:
        ada_required = 5
    elif required <= 200:
        ada_required = 6
    else:
        ada_required = max(6, int(math.ceil(required * 0.02)))
    van_accessible = int(math.ceil(ada_required / 6)) if ada_required > 0 else 0
    return required, ada_required, van_accessible

# Recompute drawn areas
building_area_drawn = 0.0
parking_lots_area = 0.0
for idx, feat in enumerate(st.session_state["workspace_features"]):
    role = st.session_state["role_map"].get(str(idx), "Unclassified")
    geom = feat.get("geometry")
    area = geodesic_area_sqft(geom)
    if role == "Building Footprint":
        building_area_drawn += area
    elif role == "Parking Lot":
        parking_lots_area += area

if building_area_drawn > 0:
    building_footprint_sf = int(round(building_area_drawn))
else:
    building_footprint_sf = int(round(building_footprint_target))

parking_est_from_area = int(math.floor(parking_lots_area / STANDARD_SPACE_SF)) if parking_lots_area > 0 else 0
drawn_parking_count = parking_est_from_area
if prefer_drawn_parking and drawn_parking_count > 0:
    provided_parking = drawn_parking_count
else:
    provided_parking = int(available_parking_input)

available_parking_area = max(0.0, site_area_sf - building_footprint_sf)

required_parking, ada_required, van_accessible = calculate_parking_requirement(studio_units, br1_units, br2_units, br3_units, project_type, state)
required_standard = max(0, required_parking - ada_required)
required_parking_area = ((required_standard * STANDARD_SPACE_SF) + max(0, (ada_required - van_accessible) * ACCESSIBLE_SPACE_SF) + (van_accessible * VAN_SPACE_SF))
total_required_area = required_parking_area * 1.3

parking_deficit = required_parking - provided_parking
area_fits = available_parking_area >= total_required_area
area_deficit = max(0.0, total_required_area - available_parking_area)

req = fire_lane_requirements[state]
fire_lane_deficit = req["min"] - site_width

# ============ QAP SCORING CALCULATION (expanded) ============
qap_points = 0
breakdown = []
max_possible_points = 0

# Tier 1: Location & Access
if transit_access:
    qap_points += 10
    breakdown.append("Transit access: +10")
max_possible_points += 10

if opportunity_zone:
    qap_points += 15
    breakdown.append("Opportunity Zone: +15")
max_possible_points += 15

if high_opportunity_area:
    qap_points += 15
    breakdown.append("High opportunity area: +15")
max_possible_points += 15

walk_points = min(5, walk_score // 20)
qap_points += walk_points
breakdown.append(f"Walkability: +{walk_points}")
max_possible_points += 5

service_points = min(10, nearby_services * 2)
qap_points += service_points
breakdown.append(f"Essential services: +{service_points}")
max_possible_points += 10

if school_proximity:
    qap_points += 5
    breakdown.append("School proximity: +5")
max_possible_points += 5

if job_centers:
    qap_points += 10
    breakdown.append("Job centers: +10")
max_possible_points += 10

if grocery_store:
    qap_points += 5
    breakdown.append("Grocery store: +5")
max_possible_points += 5

if medical_facilities:
    qap_points += 5
    breakdown.append("Medical facilities: +5")
max_possible_points += 5

if parks_recreation:
    qap_points += 3
    breakdown.append("Parks/recreation: +3")
max_possible_points += 3

# Tier 2: Site Readiness
if site_control_documented:
    qap_points += 5
    breakdown.append("Site control: +5")
max_possible_points += 5

if environmental_clear:
    qap_points += 8
    breakdown.append("Environmental clear: +8")
max_possible_points += 8

if utilities_available:
    qap_points += 3
    breakdown.append("Utilities available: +3")
max_possible_points += 3

if zoning_compliant:
    qap_points += 5
    breakdown.append("Zoning compliant: +5")
max_possible_points += 5

# Tier 3: Site Development
if green_space:
    qap_points += 5
    breakdown.append("Green space: +5")
max_possible_points += 5

if playground:
    qap_points += 3
    breakdown.append("Playground: +3")
max_possible_points += 3

if community_building:
    qap_points += 10
    breakdown.append("Community building: +10")
max_possible_points += 10

if bike_storage:
    qap_points += 2
    breakdown.append("Bike storage: +2")
max_possible_points += 2

# Tier 4: Sustainability
energy_points_map = {
    "None": 0,
    "Energy Star": 5,
    "LEED Certified": 10,
    "LEED Silver": 12,
    "LEED Gold+": 15
}
energy_pts = energy_points_map.get(energy_certification, 0)
if energy_pts > 0:
    qap_points += energy_pts
    breakdown.append(f"Energy certification: +{energy_pts}")
max_possible_points += 15

if solar_ready:
    qap_points += 8
    breakdown.append("Solar ready: +8")
max_possible_points += 8

if water_conservation:
    qap_points += 3
    breakdown.append("Water conservation: +3")
max_possible_points += 3

if ev_charging:
    qap_points += 3
    breakdown.append("EV charging: +3")
max_possible_points += 3

# Tier 5: Community
if community_support:
    qap_points += 5
    breakdown.append("Community support: +5")
max_possible_points += 5

if local_government:
    qap_points += 8
    breakdown.append("Local government support: +8")
max_possible_points += 8

if service_coordination:
    qap_points += 5
    breakdown.append("Service coordination: +5")
max_possible_points += 5

# Tier 6: Population
if special_needs_commitment:
    qap_points += 10
    breakdown.append("Special needs commitment: +10")
max_possible_points += 10

# Parking efficiency
if parking_deficit <= 0:
    qap_points += 5
    breakdown.append("Parking efficient: +5")
elif parking_deficit > 20:
    qap_points -= 5
    breakdown.append("Parking deficit: -5")
max_possible_points += 5

# FHA compliance
if fha_compliant:
    qap_points += 3
    breakdown.append("FHA compliant: +3")
max_possible_points += 3

# Normalize to 100 (avoid division by zero)
if max_possible_points <= 0:
    normalized_score = 0
else:
    normalized_score = min(100, int((qap_points / max_possible_points) * 100))

# ============ RISK ASSESSMENT (uses normalized_score) ============
risks = []

if parking_deficit > 0:
    severity = "HIGH" if parking_deficit > 20 else "MEDIUM"
    risks.append({
        'category': 'Site Compliance',
        'issue': f'Parking deficit: {parking_deficit} spaces',
        'severity': severity,
        'impact': 'QAP rejection or point deduction',
        'probability': 'High',
        'mitigation_primary': 'Request parking variance or structured parking',
        'timeline': '3-4 months (variance) / 12-18 months (structured)',
        'cost': f'$5,000 (variance) / ${parking_deficit * 25000:,.0f} (structured)'
    })

if area_deficit > 0:
    risks.append({
        'category': 'Site Constraints',
        'issue': f'Parking area deficit: {area_deficit:,.0f} SF',
        'severity': 'MEDIUM',
        'impact': 'Cannot physically fit required parking',
        'probability': 'High',
        'mitigation_primary': 'Reduce building footprint or structured parking',
        'timeline': 'Design revision / 12-18 months',
        'cost': 'Design cost / Construction premium'
    })

if fire_lane_deficit > 0:
    risks.append({
        'category': 'Site Compliance',
        'issue': f'Fire lane deficit: {fire_lane_deficit} ft',
        'severity': 'HIGH',
        'impact': 'Building permit denial',
        'probability': 'Very High',
        'mitigation_primary': 'Fire marshal variance or site redesign',
        'timeline': '2-3 months (variance)',
        'cost': '$5,000-$10,000 (variance)'
    })

if not fha_compliant:
    risks.append({
        'category': 'Federal Compliance',
        'issue': f'FHA accessibility: {fha_min_accessible - accessible_units} units short',
        'severity': 'HIGH',
        'impact': 'Fair Housing Act violation',
        'probability': 'Certain',
        'mitigation_primary': f'Increase accessible units to {fha_min_accessible}',
        'timeline': 'Immediate design adjustment',
        'cost': 'Minimal (design only)'
    })

if not site_control_documented:
    risks.append({
        'category': 'Development Readiness',
        'issue': 'Site control not documented',
        'severity': 'MEDIUM',
        'impact': 'QAP application incomplete',
        'probability': 'High',
        'mitigation_primary': 'Execute option agreement or purchase contract',
        'timeline': '2-4 weeks',
        'cost': '$5,000-$20,000'
    })

if not environmental_clear:
    risks.append({
        'category': 'Environmental',
        'issue': 'Phase I ESA not complete',
        'severity': 'MEDIUM',
        'impact': 'Financing delays',
        'probability': 'Medium',
        'mitigation_primary': 'Complete Phase I ESA',
        'timeline': '4-6 weeks',
        'cost': '$3,000-$5,000'
    })

if normalized_score < 60:
    risks.append({
        'category': 'Competitive Funding',
        'issue': f'Low QAP score: {normalized_score}/100',
        'severity': 'HIGH',
        'impact': 'Low funding probability',
        'probability': 'High',
        'mitigation_primary': 'Improve scoring factors',
        'timeline': 'Site-dependent',
        'cost': 'Variable'
    })

high_risks = [r for r in risks if r['severity'] == 'HIGH']
medium_risks = [r for r in risks if r['severity'] == 'MEDIUM']
low_risks = [r for r in risks if r['severity'] == 'LOW']

overall_risk_score = 100 - (len(high_risks) * 25 + len(medium_risks) * 10 + len(low_risks) * 5)
overall_risk_score = max(0, overall_risk_score)

# ============ MAIN UI - TABS ============
st.title("Archi Logic AI - Feasibility Analysis")
st.caption("Functional Prototype - Patent Pending")

tab1, tab2 = st.tabs([
    "Site Compliance Analysis",
    "LIHTC QAP Competitive Scoring & Risk Assessment"
])

# ============ TAB 1: SITE COMPLIANCE ============
with tab1:
    st.header("Site Compliance Analysis")
    st.subheader("Site Map & Drawing")
    
    addr_col1, addr_col2 = st.columns([4,1])
    with addr_col1:
        address = st.text_input("Address (optional)")
    with addr_col2:
        if st.button("Locate"):
            loc = geocode_address(address)
            if loc:
                st.session_state["map_center"] = (loc["lat"], loc["lon"])
                st.session_state["map_zoom"] = 18
                st.success("Located!")
            else:
                st.error("Not found")

    center_lat, center_lon = st.session_state["map_center"]
    zoom = st.session_state["map_zoom"]

    m = folium.Map(location=[center_lat, center_lon], zoom_start=zoom, control_scale=True)

    def style_by_role(feature):
        role = feature.get('properties', {}).get('role', 'Unclassified')
        if role == "Building Footprint":
            return {"color": "#1f77b4", "fillColor": "#1f77b4", "weight": 2, "fillOpacity": 0.25}
        if role == "Parking Lot":
            return {"color": "#2ca02c", "fillColor": "#2ca02c", "weight": 2, "fillOpacity": 0.25}
        return {"color": "#7f7f7f", "fillColor": "#7f7f7f", "weight": 1, "fillOpacity": 0.15}

    for idx, feat in enumerate(st.session_state["workspace_features"]):
        geom = feat.get("geometry")
        role = st.session_state["role_map"].get(str(idx), "Unclassified")
        feature_copy = {"type": "Feature", "geometry": geom, "properties": {"id": idx, "role": role}}
        folium.GeoJson(feature_copy, name=f"feature_{idx}", tooltip=f"ID {idx} — {role}", style_function=style_by_role).add_to(m)

    draw = Draw(export=True, filename='drawn.geojson', draw_options={
        'polyline': True, 'polygon': False, 'rectangle': True, 'circle': False, 'circlemarker': False, 'marker': False
    }, edit_options={'edit': True})
    draw.add_to(m)

    st.caption("Draw lines or rectangles. Lines will be converted to polygons. Click 'Add Drawn to Workspace' after drawing.")

    fret = st_folium(m, width="100%", height=map_height)

    # Drawing controls
    col_add, col_create, col_scale = st.columns([1,1,1])
    with col_add:
        if st.button("Add Drawn to Workspace"):
            new_feats = normalize_folium_return(fret)
            count_added = 0
            for nf in new_feats:
                geometry = nf.get("geometry") if isinstance(nf, dict) else None
                if not geometry:
                    continue
                if geometry.get("type") in ("LineString", "MultiLineString"):
                    polygeom = line_to_polygon_geom_with_order_corr(geometry)
                    if not polygeom:
                        continue
                    geometry = polygeom
                feat = {"type": "Feature", "geometry": geometry, "properties": {}}
                st.session_state["workspace_features"].append(feat)
                count_added += 1
            if count_added:
                st.success(f"Added {count_added} feature(s)")
            else:
                st.info("No new features found")
                
    with col_create:
        if st.button("Create exact footprint"):
            area_sf = float(building_footprint_target)
            if area_sf <= 0:
                st.error("Target must be > 0")
            else:
                area_m2 = area_sf * 0.092903
                side_m = math.sqrt(area_m2)
                lat = center_lat
                meters_per_deg_lat, meters_per_deg_lon = meters_per_degree(lat)
                delta_lat = (side_m / 2.0) / meters_per_deg_lat
                delta_lon = (side_m / 2.0) / meters_per_deg_lon
                coords = [
                    [center_lon - delta_lon, center_lat - delta_lat],
                    [center_lon + delta_lon, center_lat - delta_lat],
                    [center_lon + delta_lon, center_lat + delta_lat],
                    [center_lon - delta_lon, center_lat + delta_lat],
                    [center_lon - delta_lon, center_lat - delta_lat],
                ]
                geom = {"type": "Polygon", "coordinates": [coords]}
                feat = {"type": "Feature", "geometry": geom, "properties": {"generated": True}}
                st.session_state["workspace_features"].append(feat)
                st.success("Created exact footprint")
                
    with col_scale:
        wf_count = len(st.session_state["workspace_features"])
        if wf_count > 0:
            select_idx = st.number_input("Feature ID to scale", min_value=0, max_value=max(0, wf_count-1), value=0, step=1)
            if st.button("Scale to target"):
                try:
                    feat = st.session_state["workspace_features"][select_idx]
                    geom = feat.get("geometry")
                    new_geom = scale_geometry_to_target(geom, building_footprint_target)
                    st.session_state["workspace_features"][select_idx]["geometry"] = new_geom
                    st.success(f"Scaled feature {select_idx}")
                except Exception as e:
                    st.error(f"Failed: {e}")

    if st.button("Clear workspace"):
        st.session_state["workspace_features"] = []
        st.session_state["role_map"] = {}
        st.success("Cleared")

    # Feature classification
    st.subheader("Feature Classification")
    wf = st.session_state["workspace_features"]
    if not wf:
        st.info("Workspace empty. Draw shapes first.")
    else:
        rows = []
        for idx, feat in enumerate(wf):
            geom = feat.get("geometry")
            area_sf = geodesic_area_sqft(geom)
            role = st.session_state["role_map"].get(str(idx), "Unclassified")
            rows.append({"id": idx, "type": geom.get("type"), "area_sf": area_sf, "role": role})
        df = pd.DataFrame(rows)

        with st.form("classify_form"):
            st.write("Classify each feature as Building Footprint, Parking Lot, or Unclassified")
            for _, r in df.iterrows():
                fid = int(r["id"])
                default = st.session_state["role_map"].get(str(fid), "Unclassified")
                sel = st.selectbox(
                    f"Feature {fid} — {r['type']} — {r['area_sf']:,.0f} SF",
                    key=f"role_{fid}",
                    options=["Unclassified", "Building Footprint", "Parking Lot"],
                    index=["Unclassified", "Building Footprint", "Parking Lot"].index(default)
                )
                st.session_state["role_map"][str(fid)] = sel
            submitted = st.form_submit_button("Apply classification")
        
        if submitted:
            st.success("Classification applied")

        display_df = df.copy()
        display_df["role"] = display_df["id"].apply(lambda i: st.session_state["role_map"].get(str(i), "Unclassified"))
        st.dataframe(display_df.style.format({"area_sf": "{:,.0f}"}), height=240)

    # Compliance Results
    st.markdown("---")
    st.subheader("Compliance Analysis Results")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Units", total_units)
    with col2:
        st.metric(
            "Accessible Units",
            f"{accessible_units} ({(accessible_units/total_units*100):.1f}%)" if total_units > 0 else "0"
        )
    with col3:
        st.metric("Available Parking", provided_parking)
    with col4:
        st.metric("Site Area", f"{site_area_sf:,} SF")

    st.markdown("---")

    # Accessibility
    st.subheader("Accessibility Compliance")
    col1, col2 = st.columns(2)
    with col1:
        if fha_compliant:
            st.success("FHA COMPLIANT")
        else:
            st.error("FHA NON-COMPLIANT")
        st.write(f"Required: {fha_min_accessible} units (5% minimum)")
        st.write(f"Provided: {accessible_units} units")
        if type_a_units > 0:
            st.write(f"Type A: {type_a_units} units")
        if type_b_units > 0:
            st.write(f"Type B: {type_b_units} units")
    with col2:
        st.write("Type A Requirements:")
        st.write("• 60\" turning radius")
        st.write("• Roll-in shower (36\" × 60\")")
        st.write("• Accessible kitchen")
        st.write("• Lowered outlets")

    st.markdown("---")

    # Parking
    st.subheader("Parking Analysis")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.write("Required Parking")
        if parking_deficit > 0:
            st.error(f"{required_parking} spaces")
        else:
            st.success(f"{required_parking} spaces")
        st.write(f"ADA spaces: {ada_required}")
        st.write(f"Van accessible: {van_accessible}")
        
    with col2:
        st.write("Area Analysis")
        if area_fits:
            st.success("AREA SUFFICIENT")
        else:
            st.error("AREA INSUFFICIENT")
        st.metric("Required", f"{total_required_area:,.0f} SF")
        st.metric("Available", f"{available_parking_area:,.0f} SF")
        
    with col3:
        st.write("Compliance Status")
        if parking_deficit > 0:
            st.error(f"DEFICIT: {parking_deficit} spaces")
        else:
            st.success(f"SURPLUS: {-parking_deficit} spaces")

    st.markdown("---")

    # Fire lane
    st.subheader("Fire Lane Compliance")
    col1, col2 = st.columns(2)
    with col1:
        st.write("Requirements")
        st.write(f"Minimum Width: {req['min']} ft")
        st.write(f"Code: {req['code']}")
        st.write(f"Notes: {req['notes']}")
    with col2:
        if fire_lane_deficit <= 0:
            st.success("COMPLIANT")
            st.metric("Provided", f"{site_width} ft", delta=f"+{-fire_lane_deficit} ft buffer")
        else:
            st.error("NON-COMPLIANT")
            st.metric("Provided", f"{site_width} ft", delta=f"-{fire_lane_deficit} ft deficit", delta_color="inverse")

    st.markdown("---")

    # Overall
    st.subheader("Overall Site Compliance")
    total_issues = len([1 for r in risks if r['severity'] == 'HIGH'])

    if total_issues == 0:
        st.success("SITE COMPLIANT - All requirements satisfied")
    elif total_issues == 1:
        st.warning("CONDITIONALLY FEASIBLE - 1 critical issue")
    else:
        st.error(f"HIGH RISK - {total_issues} critical issues identified")

# ============ TAB 2: QAP SCORING + RISK ============
with tab2:
    st.header("LIHTC QAP Competitive Scoring & Risk Assessment")
    st.write(f"Site-related scoring for {state} LIHTC 9% competitive funding")
    st.markdown("---")

    col1, col2 = st.columns([1, 2])

    with col1:
        if normalized_score >= 75:
            st.success(f"{normalized_score}")
            st.success("HIGHLY COMPETITIVE")
        elif normalized_score >= 60:
            st.warning(f"{normalized_score}")
            st.warning("COMPETITIVE")
        else:
            st.error(f"{normalized_score}")
            st.error("NEEDS IMPROVEMENT")
        st.caption("Normalized score out of 100 (based on expanded QAP factors)")
        
    with col2:
        st.write("Point Breakdown (selected factors):")
        for item in breakdown:
            st.write(f"• {item}")

    st.markdown("---")

    # Score Improvement
    st.subheader("Score Improvement Opportunities")

    improvements = []

    if not transit_access:
        improvements.append({
            'action': 'Document transit access within 0.5 mile',
            'points': 10,
            'difficulty': 'Easy if available'
        })

    if not opportunity_zone:
        improvements.append({
            'action': 'Consider site in Opportunity Zone',
            'points': 15,
            'difficulty': 'Requires different site'
        })

    if walk_score < 60:
        potential = min(5, (60 // 20) - walk_points)
        if potential > 0:
            improvements.append({
                'action': 'Select more walkable location / improve pedestrian connections',
                'points': potential,
                'difficulty': 'Site-dependent'
            })

    if nearby_services < 5:
        potential = (5 - nearby_services) * 2
        if potential > 0:
            improvements.append({
                'action': f'Document {5 - nearby_services} additional essential services',
                'points': potential,
                'difficulty': 'Medium'
            })

    if parking_deficit > 0 and parking_deficit <= 10:
        improvements.append({
            'action': 'Resolve parking deficit',
            'points': 5,
            'difficulty': 'Medium'
        })

    if improvements:
        improvements.sort(key=lambda x: x['points'], reverse=True)
        for imp in improvements[:3]:
            st.info(f"{imp['action']}\n\nAdditional points: +{imp['points']}\nDifficulty: {imp['difficulty']}")
        potential_raw = qap_points + sum(i['points'] for i in improvements[:3])
        if max_possible_points > 0:
            potential_normalized = min(100, int((potential_raw / max_possible_points) * 100))
        else:
            potential_normalized = normalized_score
        st.success(f"Potential normalized score with top improvements: {potential_normalized} points")
    else:
        st.success("No obvious improvements available - score is optimized for current site")
    st.markdown("---")

    # QAP Context
    st.subheader(f"{state} QAP Context")
    defaults = lihtc_qap_defaults.get(state, {})
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Parking Reduction Allowed", "Yes" if defaults.get('parking_reduction_allowed') else "No")
    with col2:
        st.metric("Max Reduction %", f"{defaults.get('max_parking_reduction_pct')}%")
    with col3:
        st.metric("QAP Strictness", defaults.get('qap_strictness'))

    st.markdown("---")

    # Risk Dashboard
    st.subheader("Risk Summary Dashboard")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        if len(high_risks) > 0:
            st.error(f"{len(high_risks)} High Risk")
        else:
            st.success("0 High Risk")

    with col2:
        if len(medium_risks) > 0:
            st.warning(f"{len(medium_risks)} Medium Risk")
        else:
            st.success("0 Medium Risk")

    with col3:
        st.info(f"{len(low_risks)} Low Risk")

    with col4:
        if overall_risk_score >= 70:
            st.success(f"{overall_risk_score}/100")
        elif overall_risk_score >= 50:
            st.warning(f"{overall_risk_score}/100")
        else:
            st.error(f"{overall_risk_score}/100")
        st.caption("Overall Risk Score")

    st.markdown("---")

    # High Risks Detail (detailed, decision-support oriented)
    if len(high_risks) > 0:
        st.subheader("Critical Issues (Must Resolve)")

        for idx, risk in enumerate(high_risks, 1):
            with st.expander(f"Critical Issue #{idx}: {risk['issue']}", expanded=True):
                # Basic info
                col1, col2 = st.columns(2)

                with col1:
                    st.write("Problem Analysis:")
                    st.write(f"• Category: {risk['category']}")
                    st.write(f"• Impact: {risk['impact']}")
                    st.write(f"• Probability: {risk['probability']}")

                with col2:
                    st.write("Resolution Path:")
                    st.write(f"• Primary Strategy: {risk['mitigation_primary']}")
                    st.write(f"• Estimated Timeline: {risk['timeline']}")

                # Cost breakdown
                st.markdown("---")
                st.write("Cost Analysis:")

                # Parking deficit case (APPLIED UPDATED BLOCK)
                if "Parking deficit" in risk['issue']:
                    deficit_spaces = parking_deficit  # Dynamic value
                    cost_per_space = 25000  # Assumption

                    structured_total = deficit_spaces * cost_per_space

                    col_a, col_b = st.columns(2)

                    with col_a:
                        st.info(f"""
**Option 1: Parking Variance**

Request reduction approval from local jurisdiction

**Cost:** $5,000 - $10,000 (fixed)
- Application & legal fees
- Traffic study (if required)
- Independent of deficit size

**Timeline:** 3-4 months
**Success Rate:** 60-70% (varies by jurisdiction)

**Pros:**
- Low, predictable cost
- No construction needed
- Fast if approved

**Cons:**
- Approval not guaranteed
- Requires transit proximity or other justification
- May lose 5-10 QAP points
- Some jurisdictions don't allow reductions
                        """)

                    with col_b:
                        st.error(f"""
**Option 2: Structured Parking**

Build multi-level parking garage

**Cost for Your Project:** ${structured_total:,}
- **{deficit_spaces} spaces needed** (your deficit)
- × $25,000 per space (industry average)
- = ${structured_total:,} total

**Cost per space breakdown:**
- Structure (concrete/steel): $15,000
- Decking & finishes: $4,000
- Site work: $3,000
- Engineering/permits: $2,000
- Contingency (10%): $2,500
- **Total: ~$25,000/space**

**Timeline:** 12-18 months
**Success Rate:** 100% (guaranteed solution)

**Pros:**
- Guaranteed compliance
- No QAP point loss
- Adds property value

**Cons:**
- **Very expensive** (${structured_total:,})
- Long construction timeline
- Ongoing maintenance costs
- May make project financially infeasible

**Note:** If your deficit were different, cost scales proportionally.
For example:
- 10 spaces → $250,000
- 50 spaces → $1,250,000
                        """)

                    # Recommendation
                    st.markdown("---")
                    st.write("Recommended Strategy:")

                    if deficit_spaces <= 10:
                        st.success(f"""
**Primary: Pursue Variance** (deficit is small)

With only {deficit_spaces} spaces short ({(deficit_spaces/required_parking*100):.1f}% deficit), 
variance approval is realistic if you can demonstrate:
- Transit proximity (within 0.5 mile)
- Walkable neighborhood
- Shared parking potential
- Reduced parking demand (senior housing, urban location)

**Fallback:** Structured parking (${structured_total:,})
                        """)
                    elif deficit_spaces <= 20:
                        st.warning(f"""
**Primary: Pursue Variance** (worth attempting)

{deficit_spaces} spaces = {(deficit_spaces/required_parking*100):.1f}% deficit. 
Variance possible but requires strong justification:
- Transit access
- Shared parking agreement
- Alternative transportation plan

**Fallback:** Reduce unit count or structured parking

**Alternative:** Reduce total units to match available parking 
(would need to reduce by ~{math.ceil(deficit_spaces / 1.5)} units)
                        """)
                    else:
                        st.error(f"""
**Primary: Reduce Unit Count or Find Alternative Site**

{deficit_spaces} spaces = {(deficit_spaces/required_parking*100):.1f}% deficit is too large.

**Structured parking cost (${structured_total:,}) likely makes project infeasible.**

**Variance approval unlikely** with this large a deficit.

**Realistic options:**
1. Reduce unit count by ~{math.ceil(deficit_spaces / 1.5)} units
2. Find larger site
3. Consider alternative site entirely

**Structured parking should be last resort** given cost impact.
                        """)

                # Fire lane case
                elif "Fire lane" in risk['issue']:
                    col_a, col_b = st.columns(2)

                    with col_a:
                        st.info(
                            "Option 1: Fire Marshal Variance\n\n"
                            "Cost: $5,000 - $10,000\n"
                            "- Variance application\n- Fire access plan\n- Consulting fees\n\n"
                            "Timeline: 2-3 months"
                        )

                    with col_b:
                        st.warning(
                            "Option 2: Site Redesign\n\n"
                            "Cost: $15,000 - $50,000\n"
                            "- Reduce building footprint\n- Reconfigure access\n- Redesign fees\n- Lost rental income\n\n"
                            "Timeline: 1-2 months"
                        )

                # FHA case
                elif "FHA" in risk['issue'] or "FHA accessibility" in risk['issue']:
                    units_short = max(0, fha_min_accessible - accessible_units)
                    st.info(
                        f"Solution: Increase Accessible Units\n\n"
                        f"Cost: Minimal (design modification only)\n"
                        f"- Convert {units_short} standard unit(s) to accessible\n"
                        "- No structural cost if done during design phase\n"
                        "- May require unit mix adjustment\n\n"
                        "Timeline: Immediate (design phase)"
                    )

                # Generic case
                else:
                    st.write(f"Estimated cost: {risk.get('cost', 'Variable')}")

                # Why critical
                st.markdown("---")
                st.write("Why This Blocks QAP Application:")

                if "Parking" in risk['issue']:
                    shortfall_pct = (parking_deficit / required_parking * 100) if required_parking > 0 else 0
                    st.warning(
                        f"Parking deficit = {parking_deficit} spaces ({shortfall_pct:.0f}% shortfall)\n\n"
                        "QAP Impact:\n"
                        "- Automatic point deduction (typically -10 to -20 points)\n"
                        "- May trigger threshold failure if >30% deficit\n"
                        "- Requires variance documentation or structured parking commitment\n"
                        "- Lenders require parking compliance for financing approval\n\n"
                        "Bottom Line: Cannot submit a competitive LIHTC application with this deficit unresolved."
                    )

                elif "Fire lane" in risk['issue']:
                    st.warning(
                        "Building Code Requirement:\n\n"
                        "- Fire marshal approval required for building permit.\n\n"
                        "QAP Impact:\n"
                        "- QAP requires evidence of permit-ability\n"
                        "- Cannot close financing without permits\n"
                        "- Critical path item for development timeline\n\n"
                        "Bottom Line: Project cannot proceed without fire access compliance."
                    )

                elif "FHA" in risk['issue'] or "FHA accessibility" in risk['issue']:
                    st.warning(
                        "Federal Fair Housing Act:\n\n"
                        "- 5% accessible units = federal requirement (non-negotiable).\n\n"
                        "QAP Impact:\n"
                        "- Automatic application rejection if non-compliant\n"
                        "- HUD review will catch this during underwriting\n\n"
                        "Bottom Line: Must be corrected before any application."
                    )

    st.markdown("---")

    # Medium Risks (kept concise)
    if len(medium_risks) > 0:
        st.subheader("Medium Risk Items")
        for risk in medium_risks:
            with st.expander(f"{risk['issue']}"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"Category: {risk['category']}")
                    st.write(f"Impact: {risk['impact']}")
                with col2:
                    st.write(f"Mitigation: {risk['mitigation_primary']}")
                    st.write(f"Timeline: {risk['timeline']}")

    st.markdown("---")

    # Timeline & Final Feasibility
    st.subheader("Timeline to Application-Ready")
    if len(high_risks) > 0:
        st.error("CRITICAL RISKS IDENTIFIED - NOT APPLICATION-READY\nHigh-risk item(s) must be resolved before QAP submission.\nEstimated timeline: 4-6 months")
    elif len(medium_risks) > 0:
        st.warning("MEDIUM RISKS PRESENT - CONDITIONAL READINESS\nSome items should be addressed. Estimated timeline: 6-10 weeks")
    else:
        st.success("PROJECT IS APPLICATION-READY\nNo critical risks identified. Site appears suitable for LIHTC development.")

    st.markdown("---")

    st.subheader("Final Feasibility Assessment")
    total_issues = len(high_risks) + len(medium_risks)
    if total_issues == 0:
        st.success(f"SITE FEASIBLE - PROCEED WITH DEVELOPMENT\n\nQAP Score: {normalized_score}/100\nRisk Score: {overall_risk_score}/100")
    elif total_issues <= 2:
        st.warning(f"CONDITIONALLY FEASIBLE\nQAP Score: {normalized_score}/100\nRisk Score: {overall_risk_score}/100\n{total_issues} issue(s) requiring resolution.")
    else:
        st.error(f"HIGH RISK - NOT RECOMMENDED\nQAP Score: {normalized_score}/100\nRisk Score: {overall_risk_score}/100\n{total_issues} critical issues identified.")

# ============ FOOTER ============
st.markdown("---")
st.caption(
    f"Analysis Based On: {state} building codes, FHA/ADA requirements, {state} QAP criteria\n"
    f"Code References: {req['code']}, Fair Housing Act, ADA Standards\n"
    "Disclaimer: Preliminary feasibility tool. Consult local codes, QAP requirements, and qualified professionals for final determination.\n"
    f"© {COPYRIGHT_YEAR} {AUTHOR} | {CONTACT_EMAIL} | Archi Logic AI | Patent Pending"
)