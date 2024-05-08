import sys
from jinja2 import Template

# Validate the number of arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <input_vcf_path> <report_html_path>")
    sys.exit(1)

input_vcf_path = sys.argv[1]
report_html_path = sys.argv[2]

def get_vep_impact_score(vep_impacts):
    if "HIGH" in vep_impacts and "MODERATE" in vep_impacts:
        return 100
    elif "HIGH" in vep_impacts:
        return 100
    elif "MODERATE" in vep_impacts:
        return 50
    else:
        return 10

def get_clinical_significance_score(clinical_significance):
    if "pathogenic" in clinical_significance and "likely_pathogenic" in clinical_significance:
        return 100
    elif "pathogenic" in clinical_significance:
        return 100
    elif "likely_pathogenic" in clinical_significance:
        return 50
    else:
        return 10

def get_review_status_score(review_status):
    scores = {
        "reviewed_by_expert_panel": 50,
        "criteria_provided,_multiple_submitters,_no_conflicts": 30,
        "criteria_provided,_single_submitter": 10,
        "no_assertion_for_the_criteria_provided": 5
    }
    return scores.get(review_status.lower().replace(" ", "_"), 0)

data = []

try:
    with open(input_vcf_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue  # Skip header lines
            parts = line.strip().split("\t")
            info_field = parts[7]

            # Extract VEP impacts from CSQ field
            vep_impact_data = next((field.split("=")[1] for field in info_field.split(";") if field.startswith("CSQ")), "")
            vep_impacts = set([impact.split("|")[1] for impact in vep_impact_data.split(",")])

            # Extract ClinVar data
            clinical_significances = info_field.split("|")[23].split('&')
            CLNREVSTAT = next((field.split("=")[1] for field in info_field.split(";") if "CLNREVSTAT" in field), "not_provided")
            rsID = info_field.split("|")[17]
            clndn_value = next((field.split("=")[1] for field in info_field.split(";") if "CLNDN" in field), "not_provided")
            genotype = parts[9].split(":")[0]

            # Calculate scores
            vep_score = get_vep_impact_score(vep_impacts)
            clinvar_score = get_clinical_significance_score(max(clinical_significances, key=lambda x: get_clinical_significance_score(x)))
            review_status_score = get_review_status_score(CLNREVSTAT)
            total_score = vep_score + clinvar_score + review_status_score

            data.append({
                "rsID": rsID,
                "genotype": genotype,
                "clndn": clndn_value,
                "vep_impact": ", ".join(vep_impacts),
                "clinical_significance": ", ".join(clinical_significances),
                "review_status": CLNREVSTAT,
                "score": total_score
            })

    # Sort data by score in descending order
    data.sort(key=lambda x: x['score'], reverse=True)

except FileNotFoundError:
    print(f"Error: File '{input_vcf_path}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)

# Create an HTML report
template = Template("""
<!DOCTYPE html>
<html>
<head>
    <title>Clinical Report</title>
    <style>
        table {
            border-collapse: collapse;
            width: 80%;
            margin: 20px auto;
        }
        table, th, td {
            border: 1px solid black;
        }
        th, td {
            padding: 8px;
            text-align: left;
        }
        th {
            background-color: #f2f2f2;
        }
    </style>
</head>
<body>
    <h1>🧬 Clinical Report</h1>

    <label for="hpoInput">Filter by HPO Term:</label>
    <input type="text" id="hpoInput">
    <button onclick="filterVariants()">Filter</button>
    
    <label for="reviewStatusInput">Filter by Review Status:</label>
    <input type="text" id="reviewStatusInput">
    <button onclick="filterByReviewStatus()">Filter</button>

    <table border="1" id="variantsTable">
        <tr>
            <th>Score</th>
            <th>VEP Impact</th>
            <th>Clinical Significance</th>
            <th>rsID</th>
            <th>Genotype</th>
            <th>ClinVar Annotation (CLNDN)</th>
            <th>Review Status</th>
        </tr>
        {% for variant in data %}
        <tr>
            <td>{{ variant.score }}</td>
            <td>{{ variant.vep_impact }}</td>
            <td>{{ variant.clinical_significance }}</td>
            <td>{{ variant.rsID }}</td>
            <td>{{ variant.genotype }}</td>
            <td>{{ variant.clndn }}</td>
            <td>{{ variant.review_status }}</td>
        </tr>
        {% endfor %}
    </table>

    <script>
        var variants = {{ data|tojson|safe }};
        
        function displayVariants() {
            var table = document.getElementById("variantsTable");
            for (var i = 0; i < variants.length; i++) {
                var variant = variants[i];
                var row = table.insertRow(-1);
                row.insertCell(0).innerHTML = variant.score;
                row.insertCell(1).innerHTML = variant.vep_impact;
                row.insertCell(2).innerHTML = variant.clinical_significance;
                row.insertCell(3).innerHTML = variant.rsID;
                row.insertCell(4).innerHTML = variant.genotype;
                row.insertCell(5).innerHTML = variant.clndn;
                row.insertCell(6).innerHTML = variant.review_status || "Not available";
            }
        }

        function filterVariants() {
            var hpoFilter = document.getElementById("hpoInput").value.toLowerCase();
            var table = document.getElementById("variantsTable");
            var rows = table.getElementsByTagName("tr");

            for (var i = 1; i < rows.length; i++) { // Start from 1 to skip header row
                var row = rows[i];
                var clndn = row.cells[5].textContent.toLowerCase();

                if (clndn.includes(hpoFilter)) {
                    row.style.display = "";
                } else {
                    row.style.display = "none";
                }
            }
        }

        function filterByReviewStatus() {
            var reviewStatusFilter = document.getElementById("reviewStatusInput").value.toLowerCase();
            var table = document.getElementById("variantsTable");
            var rows = table.getElementsByTagName("tr");

            for (var i = 1; i < rows.length; i++) { // Start from 1 to skip header row
                var row = rows[i];
                var reviewStatus = row.cells[6].textContent.toLowerCase();

                if (reviewStatus.includes(reviewStatusFilter)) {
                    row.style.display = "";
                } else {
                    row.style.display = "none";
                }
            }
        }

        // Initialize the report by displaying all variants
        displayVariants();
    </script>
</body>
</html>
""")

# Write the HTML report to the output file
try:
    with open(report_html_path, "w") as report_file:
        report_file.write(template.render(data=data))
except Exception as e:
    print(f"An error occurred while writing the report: {e}")
